#!/usr/bin/env python3
"""
04_multicopy_coding/analyze.py
Task: Multi-copy genes with high homology in coding regions.

Context: Some gene families have multiple functional copies (paralogs) with
highly similar CDS sequences.  Classic examples: HLA-A/B/C, olfactory receptors,
immunoglobulin constant genes, ribosomal protein paralogs.  For these genes the
CDS alone is insufficient to achieve unique mapping — the 3' UTR may provide the
only distinguishing sequence.

This script:
  1. Identifies protein-coding genes in 'local_tandem_dup' or 'distant_paralog'
     phenomena from the gene overlap annotation — these are sequence-cluster
     members where the cluster contains multiple PC genes.
  2. For MANE-level genes in those clusters, extracts STAR UF (overall vs 3'
     window) to quantify whether the UTR rescues uniqueness.
  3. Computes cluster-level statistics: n members, mean UF, frac rescued by UTR.
  4. Flags clusters where even the UTR provides no rescue (worst case for scRNA-seq).

Ambiguity type: *Paralog ambiguity* (conserved CDS across gene family members)

Inputs (independent):
  ../sc3UTRs/results/gene_overlap_annotation.tsv
  ../sc3UTRs/results/sc3prime_v5g_MANE_RNA_PE_100.tsv
  ../../Ref/MANE.GRCh38.v1.5.ensembl_rna.fna
  ../results/id_conversion_table.tsv

Outputs:
  multicopy_genes.tsv        per-gene UF and UTR rescue status
  multicopy_clusters.tsv     per-sequence-cluster: n_members, UF, UTR rescue fraction
  utr_rescue_summary.tsv     aggregate rescue rates by n_cluster_members
"""

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

HERE   = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from lib_check import check_library
ROOT   = HERE.parent.parent.parent
REFDIR = ROOT / 'Ref'
MAPDIR = ROOT / 'mappability'

ANNOT      = MAPDIR / 'sc3UTRs' / 'results' / 'gene_overlap_annotation.tsv'
SC3P_MANE  = MAPDIR / 'sc3UTRs' / 'results' / 'sc3prime_v5g_MANE_RNA_PE_100.tsv'
MANE_FASTA = REFDIR / 'MANE.GRCh38.v1.5.ensembl_rna.fna'
ID_TABLE   = MAPDIR / 'results' / 'id_conversion_table.tsv'

# Phenomena where sequence-cluster membership drives the issue
TARGET_PHENOMENA = {'local_tandem_dup', 'distant_paralog'}
# Also include protein-coding genes in clusters with seq_cluster_n_pc > 1
# regardless of primary phenomenon label


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


def utr_rescue(uf_overall, uf_3p250, threshold=0.1):
    """
    'strong'  — uf_3p250 > uf_overall + threshold  (UTR adds clear uniqueness)
    'partial' — uf_3p250 > uf_overall               (small improvement)
    'none'    — uf_3p250 <= uf_overall              (no UTR help)
    'lost'    — uf_3p250 < uf_overall - threshold   (3' UTR overlaps something else)
    """
    d = uf_3p250 - uf_overall
    if d > threshold:
        return 'strong'
    elif d > 0:
        return 'partial'
    elif d > -threshold:
        return 'none'
    else:
        return 'lost'


def main():
    # ── Library consistency check ─────────────────────────────────────────────
    # Single input run: sc3prime_v5g_MANE_RNA_PE_100.tsv (MANE RNA PE 100 bp).
    # UF values are fractions [0,1]; no cross-run mixing here.
    check_library(
        run_label   = 'v5g_MANE_RNA_PE_100',
        read_len_bp = 100,
        mode        = 'PE',
    )

    print('Loading gene overlap annotation …')
    annot = load_annotation(ANNOT)

    print('Loading MANE FASTA gene map …')
    fasta_map = load_mane_fasta_gene_map(MANE_FASTA)

    print('Loading STAR UF — MANE RNA PE 100 bp …')
    sc3p = load_sc3prime(SC3P_MANE)

    # ── Identify multicopy genes ──────────────────────────────────────────────
    # A gene qualifies if:
    #   - phenomenon in TARGET_PHENOMENA (local_tandem_dup, distant_paralog)
    # OR
    #   - protein-coding, seq_cluster_n_pc > 1  (in a PC-only or mixed cluster)
    multicopy_gids = {}
    for gid, row in annot.items():
        if row['gene_type'] != 'protein_coding':
            continue
        try:
            n_pc_in_cluster = int(row.get('seq_cluster_n_pc', 0))
        except ValueError:
            n_pc_in_cluster = 0
        if row['phenomenon'] in TARGET_PHENOMENA or n_pc_in_cluster > 1:
            multicopy_gids[gid] = row

    print(f'  {len(multicopy_gids):,} protein-coding multi-copy genes')

    # Build gene_symbol → annot for MANE lookup
    sym2annot = defaultdict(list)
    for gid, row in multicopy_gids.items():
        sym2annot[row['gene_name']].append(row)

    # ── Build per-gene output rows ────────────────────────────────────────────
    gene_rows = []
    for enst, gsym in sorted(fasta_map.items(), key=lambda x: x[1]):
        arow_list = sym2annot.get(gsym, [])
        if not arow_list:
            continue
        arow = arow_list[0]   # take first (protein_coding preferred above)

        sc3row = sc3p.get(enst, {})
        ufo    = sc3row.get('uf_overall', '')
        uf3    = sc3row.get('uf_3p250', '')
        dlt    = sc3row.get('delta_3p250', '')

        rescue = ''
        if ufo != '' and uf3 != '':
            rescue = utr_rescue(float(ufo), float(uf3))

        gene_rows.append({
            'enst_id':              enst,
            'gene_symbol':          gsym,
            'phenomenon':           arow['phenomenon'],
            'seq_cluster_id':       arow.get('seq_cluster_id', ''),
            'seq_cluster_size':     arow.get('seq_cluster_size', ''),
            'seq_cluster_n_genes':  arow.get('seq_cluster_n_genes', ''),
            'seq_cluster_n_pc':     arow.get('seq_cluster_n_pc', ''),
            'seq_cluster_mixed':    arow.get('seq_cluster_mixed', ''),
            'chrom':                arow.get('chrom', ''),
            'seq_cluster_max_span_mb': arow.get('seq_cluster_max_span_mb', ''),
            'uf_overall':           ufo,
            'uf_3p250':             uf3,
            'delta_3p250':          dlt,
            'utr_rescue':           rescue,
        })

    gene_fields = ['enst_id', 'gene_symbol', 'phenomenon',
                   'seq_cluster_id', 'seq_cluster_size', 'seq_cluster_n_genes',
                   'seq_cluster_n_pc', 'seq_cluster_mixed', 'chrom',
                   'seq_cluster_max_span_mb',
                   'uf_overall', 'uf_3p250', 'delta_3p250', 'utr_rescue']
    with open(HERE / 'multicopy_genes.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=gene_fields, delimiter='\t', lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(gene_rows)
    print(f'\nSaved multicopy_genes.tsv  ({len(gene_rows)} genes)')

    # ── Cluster-level summary ─────────────────────────────────────────────────
    clusters = defaultdict(list)
    for r in gene_rows:
        cid = r['seq_cluster_id']
        if cid:
            clusters[cid].append(r)

    clust_rows = []
    for cid, members in sorted(clusters.items()):
        ufs  = [float(m['uf_overall']) for m in members if m['uf_overall'] != '']
        uf3s = [float(m['uf_3p250'])   for m in members if m['uf_3p250'] != '']
        rescues = [m['utr_rescue'] for m in members if m['utr_rescue']]
        clust_rows.append({
            'seq_cluster_id':       cid,
            'n_members_mane':       len(members),
            'seq_cluster_n_pc':     members[0].get('seq_cluster_n_pc', ''),
            'seq_cluster_n_genes':  members[0].get('seq_cluster_n_genes', ''),
            'seq_cluster_mixed':    members[0].get('seq_cluster_mixed', ''),
            'span_mb':              members[0].get('seq_cluster_max_span_mb', ''),
            'gene_symbols':         ';'.join(m['gene_symbol'] for m in members),
            'chroms':               ';'.join(sorted({m['chrom'] for m in members})),
            'mean_uf_overall':      round(mean(ufs), 4) if ufs else '',
            'min_uf_overall':       round(min(ufs), 4) if ufs else '',
            'frac_uf_zero':         round(sum(1 for u in ufs if u == 0) / len(ufs), 4) if ufs else '',
            'mean_uf_3p250':        round(mean(uf3s), 4) if uf3s else '',
            'min_uf_3p250':         round(min(uf3s), 4) if uf3s else '',
            'n_strong_rescue':      rescues.count('strong'),
            'n_partial_rescue':     rescues.count('partial'),
            'n_no_rescue':          rescues.count('none'),
            'n_lost_in_utr':        rescues.count('lost'),
            'worst_case':           (
                'utr_also_lost'  if rescues.count('lost') > 0 else
                'no_utr_rescue'  if rescues.count('none') == len(rescues) and rescues else
                'partial_rescue' if rescues.count('strong') == 0 and rescues else
                'rescued'
            ),
        })

    clust_fields = ['seq_cluster_id', 'n_members_mane', 'seq_cluster_n_pc',
                    'seq_cluster_n_genes', 'seq_cluster_mixed', 'span_mb',
                    'gene_symbols', 'chroms',
                    'mean_uf_overall', 'min_uf_overall', 'frac_uf_zero',
                    'mean_uf_3p250', 'min_uf_3p250',
                    'n_strong_rescue', 'n_partial_rescue', 'n_no_rescue', 'n_lost_in_utr',
                    'worst_case']
    with open(HERE / 'multicopy_clusters.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=clust_fields, delimiter='\t', lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(clust_rows)
    print(f'Saved multicopy_clusters.tsv  ({len(clust_rows)} clusters)')

    # ── Rescue summary by cluster size bin ────────────────────────────────────
    size_bins = defaultdict(lambda: {'strong': 0, 'partial': 0, 'none': 0, 'lost': 0, 'n': 0})
    for r in gene_rows:
        try:
            n_pc = int(r.get('seq_cluster_n_pc', 0))
        except ValueError:
            n_pc = 0
        n_bin = min(n_pc, 10)
        size_bins[n_bin]['n'] += 1
        resc = r.get('utr_rescue', '')
        if resc in ('strong', 'partial', 'none', 'lost'):
            size_bins[n_bin][resc] += 1

    rescue_rows = []
    for n in sorted(size_bins):
        d  = size_bins[n]
        nn = d['n']
        rescue_rows.append({
            'n_pc_in_cluster': f'{n}+' if n == 10 else str(n),
            'n_genes':         nn,
            'n_strong_rescue': d['strong'],
            'n_partial_rescue': d['partial'],
            'n_no_rescue':     d['none'],
            'n_lost_in_utr':   d['lost'],
            'frac_rescued':    round((d['strong'] + d['partial']) / nn, 4) if nn else '',
            'frac_no_rescue':  round(d['none'] / nn, 4) if nn else '',
        })

    rescue_fields = ['n_pc_in_cluster', 'n_genes', 'n_strong_rescue', 'n_partial_rescue',
                     'n_no_rescue', 'n_lost_in_utr', 'frac_rescued', 'frac_no_rescue']
    with open(HERE / 'utr_rescue_summary.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=rescue_fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rescue_rows)
    print(f'Saved utr_rescue_summary.tsv')

    # ── Console ───────────────────────────────────────────────────────────────
    all_resc = [r['utr_rescue'] for r in gene_rows if r['utr_rescue']]
    print('\n═══ Summary ═══')
    print(f'Multi-copy PC genes (MANE) : {len(gene_rows):,}')
    print(f'Sequence clusters          : {len(clust_rows):,}')
    if all_resc:
        for cat in ('strong', 'partial', 'none', 'lost'):
            n = all_resc.count(cat)
            print(f'  UTR rescue {cat:<8}    : {n:,}  ({100*n/len(all_resc):.1f}%)')
    worst = [r for r in clust_rows if r['worst_case'] in ('utr_also_lost', 'no_utr_rescue')]
    print(f'Clusters with no UTR rescue: {len(worst):,}')
    print('  Examples:', '; '.join(r['gene_symbols'][:40] for r in worst[:5]))


if __name__ == '__main__':
    main()
