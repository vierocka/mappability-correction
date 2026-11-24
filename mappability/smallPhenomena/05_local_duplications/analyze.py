#!/usr/bin/env python3
"""
05_local_duplications/analyze.py
Task: Local gene duplications — tandem or proximal copies on the same chromosome.

Context: Local duplication produces gene pairs (or clusters) that often share
not just CDS but also 5'/3' flanking sequence.  Because they arose by duplication,
their sequences are highly similar throughout.  The flanking/UTR sequence may still
be distinct enough to rescue unique mapping.  This is in contrast to distant
paralogs, where the diverged UTR usually helps.

Focus here: tandem or local duplications (seq_cluster_max_span_mb ≤ 5 Mb, same
chromosome), captured by the 'local_tandem_dup' phenomenon label.  We compare
kallisto (EM redistribution) vs STAR (discard-by-default) behaviour.

Inputs (independent):
  ../sc3UTRs/results/gene_overlap_annotation.tsv
  ../sc3UTRs/results/sc3prime_v5g_MANE_RNA_PE_100.tsv        (STAR, MANE RNA PE)
  ../kallisto_results/mappability_MANE_RNA_PE/kallisto_vs_MANE/abundance.tsv
  ../../Ref/MANE.GRCh38.v1.5.ensembl_rna.fna

Outputs:
  local_dup_genes.tsv          per-gene: UF, TPM, rescue status
  local_dup_clusters.tsv       per-cluster: STAR UF vs kallisto TPM distribution
  aligner_comparison.tsv       for each local-dup gene: star_uf vs kallisto_tpm
                                to show how the two aligners handle reads differently
"""

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

HERE   = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from lib_check import warn_if_mixing_runs
ROOT   = HERE.parent.parent.parent
REFDIR = ROOT / 'Ref'
MAPDIR = ROOT / 'mappability'

ANNOT      = MAPDIR / 'sc3UTRs' / 'results' / 'gene_overlap_annotation.tsv'
SC3P_MANE  = MAPDIR / 'sc3UTRs' / 'results' / 'sc3prime_v5g_MANE_RNA_PE_100.tsv'
KAL_ABUND  = (MAPDIR / 'kallisto_results'
              / 'mappability_MANE_RNA_PE' / 'kallisto_vs_MANE' / 'abundance.tsv')
MANE_FASTA = REFDIR / 'MANE.GRCh38.v1.5.ensembl_rna.fna'

LOCAL_SPAN_MB = 5.0     # max span (Mb) to be considered "local"


def load_annotation(path):
    rows = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            rows[row['gene_id']] = row
    return rows


def load_sc3prime(path):
    rows = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            enst = row['transcript_id'].split('.')[0]
            rows[enst] = row
    return rows


def load_kallisto(path):
    kal = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            enst = row['target_id'].split('.')[0]
            kal[enst] = float(row['tpm'])
    return kal


def load_mane_fasta_gene_map(path):
    result = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            enst  = line[1:].split()[0].split('.')[0]
            gsym  = re.search(r'gene_symbol:(\S+)', line)
            result[enst] = gsym.group(1) if gsym else ''
    return result


def main():
    # ── Library consistency check ─────────────────────────────────────────────
    # STAR UF: v5g MANE RNA PE 100 bp
    # kallisto TPM: mappability_MANE_RNA_PE — also PE 100 bp
    # Same read length AND same reference (MANE) → STAR vs kallisto comparison valid.
    # UF is a fraction; TPM is per-million.
    # IMPORTANT: the comparison here is ALIGNER BEHAVIOUR, not expression level.
    # STAR discards multi-mappers (low UF = reads lost).
    # kallisto redistributes via EM (all reads assigned, may inflate TPM).
    warn_if_mixing_runs([
        (SC3P_MANE, 100),   # STAR UF — MANE RNA PE 100 bp
        (KAL_ABUND, 100),   # kallisto TPM — MANE RNA PE 100 bp
    ])

    print('Loading gene overlap annotation …')
    annot = load_annotation(ANNOT)

    print('Loading MANE FASTA gene map …')
    fasta_map = load_mane_fasta_gene_map(MANE_FASTA)

    print('Loading STAR UF …')
    sc3p = load_sc3prime(SC3P_MANE)

    print('Loading kallisto TPM …')
    kal  = load_kallisto(KAL_ABUND)

    # ── Identify local tandem duplications ────────────────────────────────────
    # Primary criterion: phenomenon == 'local_tandem_dup'
    # Additional: seq_cluster_max_span_mb ≤ LOCAL_SPAN_MB and same-chrom cluster
    local_dup_gids = {}
    for gid, row in annot.items():
        if row['gene_type'] != 'protein_coding':
            continue
        if row['phenomenon'] == 'local_tandem_dup':
            local_dup_gids[gid] = row
        elif row['phenomenon'] == 'distant_paralog':
            # Check if actually local by span
            try:
                span = float(row.get('seq_cluster_max_span_mb', 999))
            except ValueError:
                span = 999
            multi_chrom = row.get('seq_cluster_multi_chrom', 'True').strip().lower()
            if span <= LOCAL_SPAN_MB and multi_chrom in ('false', '0', 'no'):
                local_dup_gids[gid] = row

    print(f'  {len(local_dup_gids):,} local duplication genes')

    sym2annot = defaultdict(list)
    for gid, row in local_dup_gids.items():
        sym2annot[row['gene_name']].append(row)

    # ── Build per-gene output ─────────────────────────────────────────────────
    gene_rows = []
    for enst, gsym in sorted(fasta_map.items(), key=lambda x: x[1]):
        arow_list = sym2annot.get(gsym, [])
        if not arow_list:
            continue
        arow = arow_list[0]

        sc3row  = sc3p.get(enst, {})
        kal_tpm = kal.get(enst, None)

        ufo  = sc3row.get('uf_overall', '')
        uf3  = sc3row.get('uf_3p250', '')
        dlt  = sc3row.get('delta_3p250', '')

        # STAR vs kallisto: STAR UF = fraction of reads that map uniquely (low → multi-mapping)
        # kallisto TPM = EM-redistributed expression estimate (potentially overestimated for
        # non-unique genes because EM spreads counts proportionally to prior)
        star_unique = ufo
        kal_est     = round(kal_tpm, 4) if kal_tpm is not None else ''

        gene_rows.append({
            'enst_id':               enst,
            'gene_symbol':           gsym,
            'phenomenon':            arow['phenomenon'],
            'chrom':                 arow.get('chrom', ''),
            'seq_cluster_id':        arow.get('seq_cluster_id', ''),
            'seq_cluster_n_pc':      arow.get('seq_cluster_n_pc', ''),
            'seq_cluster_max_span_mb': arow.get('seq_cluster_max_span_mb', ''),
            'star_uf_overall':       star_unique,
            'star_uf_3p250':         uf3,
            'star_delta_3p250':      dlt,
            'kallisto_tpm':          kal_est,
            'aligner_concordance': (
                # High STAR UF + high kallisto TPM → both agree, good signal
                # Low STAR UF + high kallisto TPM → kallisto inflates (EM redistribution)
                # Low STAR UF + low kallisto TPM  → gene genuinely lowly expressed or lost
                'inflation_risk' if (
                    ufo != '' and kal_est != '' and
                    float(ufo) < 0.3 and float(kal_est) > 0
                ) else
                'discarded_by_star' if (
                    ufo != '' and ufo != 'NA' and float(ufo) == 0
                ) else
                'concordant'
            ),
        })

    gene_fields = ['enst_id', 'gene_symbol', 'phenomenon', 'chrom',
                   'seq_cluster_id', 'seq_cluster_n_pc', 'seq_cluster_max_span_mb',
                   'star_uf_overall', 'star_uf_3p250', 'star_delta_3p250',
                   'kallisto_tpm', 'aligner_concordance']
    with open(HERE / 'local_dup_genes.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=gene_fields, delimiter='\t', lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(gene_rows)
    print(f'\nSaved local_dup_genes.tsv  ({len(gene_rows)} genes)')

    # ── Cluster-level: STAR vs kallisto ───────────────────────────────────────
    clusters = defaultdict(list)
    for r in gene_rows:
        cid = r['seq_cluster_id']
        if cid:
            clusters[cid].append(r)

    clust_rows = []
    for cid, members in sorted(clusters.items()):
        star_ufs = [float(m['star_uf_overall']) for m in members if m['star_uf_overall'] not in ('', 'NA')]
        uf3s     = [float(m['star_uf_3p250'])   for m in members if m['star_uf_3p250']   not in ('', 'NA')]
        tpms     = [float(m['kallisto_tpm'])     for m in members if m['kallisto_tpm']    not in ('', 'NA')]
        concords = [m['aligner_concordance'] for m in members]
        clust_rows.append({
            'seq_cluster_id':       cid,
            'n_members_mane':       len(members),
            'seq_cluster_n_pc':     members[0].get('seq_cluster_n_pc', ''),
            'span_mb':              members[0].get('seq_cluster_max_span_mb', ''),
            'chroms':               ';'.join(sorted({m['chrom'] for m in members})),
            'gene_symbols':         ';'.join(m['gene_symbol'] for m in members),
            'mean_star_uf':         round(mean(star_ufs), 4) if star_ufs else '',
            'min_star_uf':          round(min(star_ufs), 4) if star_ufs else '',
            'frac_star_uf_zero':    round(sum(1 for u in star_ufs if u == 0) / len(star_ufs), 4) if star_ufs else '',
            'mean_star_uf_3p250':   round(mean(uf3s), 4) if uf3s else '',
            'total_kallisto_tpm':   round(sum(tpms), 3) if tpms else '',
            'max_tpm_fraction':     round(max(tpms) / sum(tpms), 4) if tpms and sum(tpms) > 0 else '',
            'n_inflation_risk':     concords.count('inflation_risk'),
            'n_discarded':          concords.count('discarded_by_star'),
        })

    clust_fields = ['seq_cluster_id', 'n_members_mane', 'seq_cluster_n_pc', 'span_mb',
                    'chroms', 'gene_symbols',
                    'mean_star_uf', 'min_star_uf', 'frac_star_uf_zero',
                    'mean_star_uf_3p250', 'total_kallisto_tpm', 'max_tpm_fraction',
                    'n_inflation_risk', 'n_discarded']
    with open(HERE / 'local_dup_clusters.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=clust_fields, delimiter='\t', lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(clust_rows)
    print(f'Saved local_dup_clusters.tsv  ({len(clust_rows)} clusters)')

    # ── Aligner comparison ────────────────────────────────────────────────────
    # Every gene: star UF vs kallisto TPM  → shows discordance clearly
    aln_fields = ['enst_id', 'gene_symbol', 'seq_cluster_id', 'seq_cluster_n_pc',
                  'star_uf_overall', 'star_uf_3p250', 'kallisto_tpm', 'aligner_concordance']
    with open(HERE / 'aligner_comparison.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=aln_fields, delimiter='\t', lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader()
        w.writerows(sorted(gene_rows, key=lambda r: r['gene_symbol']))
    print(f'Saved aligner_comparison.tsv')

    # ── Console ───────────────────────────────────────────────────────────────
    concords = [r['aligner_concordance'] for r in gene_rows]
    star_ufs = [float(r['star_uf_overall']) for r in gene_rows if r['star_uf_overall'] not in ('', 'NA')]
    print('\n═══ Summary ═══')
    print(f'Local dup genes (MANE)  : {len(gene_rows):,}')
    print(f'Clusters                : {len(clust_rows):,}')
    if star_ufs:
        print(f'Mean STAR UF            : {mean(star_ufs):.4f}')
        print(f'Frac STAR UF = 0        : {sum(1 for u in star_ufs if u==0)/len(star_ufs):.3f}')
    for cat in ('concordant', 'inflation_risk', 'discarded_by_star'):
        n = concords.count(cat)
        print(f'  {cat:<25}: {n:,}  ({100*n/max(len(concords),1):.1f}%)')


if __name__ == '__main__':
    main()
