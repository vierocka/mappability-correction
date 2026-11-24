#!/usr/bin/env python3
"""
02_overlapping_annotation/analyze.py
Task: Genes with overlapping genomic annotation at MANE level.

Context: GENCODE v49 annotates ~78k genes.  Many share genomic space: one gene's
exons overlap another gene on the same or opposite strand.  This directly causes
multi-mapping in RNA-seq because a read spanning the overlap zone is compatible
with both transcripts.

This script:
  1. Loads the MANE transcript set (GFF + RNA FASTA)
  2. Loads overlap classification from gene_overlap_annotation.tsv
  3. Loads STAR UF from sc3prime_v5g_MANE_RNA_PE_100.tsv  (canonical: MANE RNA PE 100 bp)
  4. Computes per-gene overlap metrics: what fraction of the transcript is in
     the overlap zone (estimated as seq_cluster or overlap group membership)
  5. Pulls kallisto TPM from mappability_MANE_RNA_PE/kallisto_vs_MANE/abundance.tsv
  6. Reports: how are reads assigned across overlapping gene pairs?

Ambiguity type: *Paralog ambiguity* (same-strand overlaps, e.g. tandem genes)
                *Genomic overlap* (antisense overlaps, nested genes)

Inputs (independent — no dependency on other smallPhenomena tasks):
  ../../Ref/MANE.GRCh38.v1.5.ensembl_rna.fna
  ../../Ref/MANE.GRCh38.v1.5.ensembl_genomic.gff
  ../sc3UTRs/results/gene_overlap_annotation.tsv
  ../sc3UTRs/results/sc3prime_v5g_MANE_RNA_PE_100.tsv
  ../kallisto_results/mappability_MANE_RNA_PE/kallisto_vs_MANE/abundance.tsv
  ../results/id_conversion_table.tsv

Outputs:
  mane_overlap_genes.tsv        MANE genes with overlap phenomena + UF + kallisto TPM
  overlap_group_summary.tsv     Per overlap group: genes, UF stats, read distribution
  phenomenon_uf_summary.tsv     Mean UF per phenomenon type (antisense/sense/nested/unique)
"""

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

HERE   = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from lib_check import check_library, warn_if_mixing_runs
ROOT   = HERE.parent.parent.parent
REFDIR = ROOT / 'Ref'
MAPDIR = ROOT / 'mappability'

MANE_FASTA = REFDIR / 'MANE.GRCh38.v1.5.ensembl_rna.fna'
MANE_GFF   = REFDIR / 'MANE.GRCh38.v1.5.ensembl_genomic.gff'
ANNOT      = MAPDIR / 'sc3UTRs' / 'results' / 'gene_overlap_annotation.tsv'
SC3PRIME   = MAPDIR / 'sc3UTRs' / 'results' / 'sc3prime_v5g_MANE_RNA_PE_100.tsv'
KAL_ABUND  = (MAPDIR / 'kallisto_results'
              / 'mappability_MANE_RNA_PE' / 'kallisto_vs_MANE' / 'abundance.tsv')
ID_TABLE   = MAPDIR / 'results' / 'id_conversion_table.tsv'


# ── loaders ───────────────────────────────────────────────────────────────────

def load_mane_fasta_lengths(path):
    """Returns {ENST_no_version: (gene_symbol, transcript_len)}."""
    result = {}
    cur_id = None
    cur_len = 0
    cur_sym = ''
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if cur_id:
                    result[cur_id] = (cur_sym, cur_len)
                parts = line.split()
                enst  = parts[0][1:].split('.')[0]
                gsym_m = re.search(r'gene_symbol:(\S+)', line)
                cur_id  = enst
                cur_sym = gsym_m.group(1) if gsym_m else ''
                cur_len = 0
            else:
                cur_len += len(line.strip())
    if cur_id:
        result[cur_id] = (cur_sym, cur_len)
    return result


def load_id_table(path):
    """Returns {ENST_no_version: gene_symbol}."""
    result = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            enst = row.get('ensembl_transcript_id', '').split('.')[0]
            if enst:
                result[enst] = row.get('gene_symbol', '')
    return result


def load_annotation(path):
    """Returns {gene_id: row_dict} from gene_overlap_annotation.tsv."""
    rows = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            rows[row['gene_id']] = row
    return rows


def load_sc3prime(path):
    """Returns {ENST_no_version: row_dict}."""
    rows = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            enst = row['transcript_id'].split('.')[0]
            rows[enst] = row
        return rows


def load_kallisto(path):
    """Returns {ENST_no_version: tpm}."""
    kal = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            enst = row['target_id'].split('.')[0]
            kal[enst] = float(row['tpm'])
    return kal


def load_mane_gff(path):
    """
    Returns {ENSG_no_version: (gene_name, chrom, start, end, strand)}
    from gene-level features only.
    """
    genes = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip().split('\t')
            if len(cols) < 9 or cols[2] != 'gene':
                continue
            attrs = cols[8]
            gid_m = re.search(r'gene_id=([^;]+)', attrs)
            gnm_m = re.search(r'gene_name=([^;]+)', attrs)
            if not gid_m:
                continue
            gid = gid_m.group(1).split('.')[0]
            genes[gid] = (
                gnm_m.group(1) if gnm_m else '',
                cols[0], int(cols[3]), int(cols[4]), cols[6]
            )
    return genes


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    # ── Library consistency check ─────────────────────────────────────────────
    # sc3prime_v5g:  MANE RNA PE 100 bp — single canonical run
    # kallisto:      MANE RNA PE 100 bp (mappability_MANE_RNA_PE)
    # Both 100 bp PE → direct comparison valid.  UF and TPM are %-normalised.
    warn_if_mixing_runs([
        (SC3PRIME,  100),   # v5g_MANE_RNA_PE_100 — STAR UF
        (KAL_ABUND, 100),   # kallisto TPM — MANE RNA PE 100 bp
    ])

    print('Loading MANE FASTA …')
    mane_lengths = load_mane_fasta_lengths(MANE_FASTA)
    print(f'  {len(mane_lengths):,} MANE transcripts')

    print('Loading ID conversion table …')
    id_table = load_id_table(ID_TABLE)

    print('Loading gene overlap annotation …')
    annot = load_annotation(ANNOT)

    print('Loading STAR UF (MANE RNA PE 100 bp) …')
    sc3p = load_sc3prime(SC3PRIME)

    print('Loading kallisto TPM (MANE RNA PE → MANE index) …')
    kal  = load_kallisto(KAL_ABUND)

    print('Loading MANE GFF gene coordinates …')
    mane_genes = load_mane_gff(MANE_GFF)

    # ── Build per-MANE-gene row ───────────────────────────────────────────────
    # Key: gene may appear in annot by ENSG, in id_table by ENST
    # Strategy: iterate over MANE transcripts (from fasta), resolve gene_id via
    # gene_name → annot lookup, then pull UF from sc3prime by ENST

    # Build gene_symbol → annotation rows (multiple genes can have same symbol)
    sym2annot = defaultdict(list)
    for gid, arow in annot.items():
        sym2annot[arow['gene_name']].append(arow)

    gene_rows = []
    for enst, (sym, tx_len) in sorted(mane_lengths.items()):
        sc3row = sc3p.get(enst, {})
        kal_tpm = kal.get(enst, None)

        # Resolve annotation: match by gene_name
        arow_candidates = sym2annot.get(sym, [])
        # Prefer protein_coding MANE gene
        arow = None
        for c in arow_candidates:
            if c.get('gene_type') == 'protein_coding':
                arow = c
                break
        if arow is None and arow_candidates:
            arow = arow_candidates[0]

        uf_overall = sc3row.get('uf_overall', '')
        uf_3p250   = sc3row.get('uf_3p250', '')
        delta_3p250 = sc3row.get('delta_3p250', '')

        row = {
            'enst_id':             enst,
            'gene_symbol':         sym,
            'transcript_len_nt':   tx_len,
            'phenomenon':          arow['phenomenon'] if arow else 'unknown',
            'overlap_group_id':    arow.get('overlap_group_id', '') if arow else '',
            'overlap_group_size':  arow.get('overlap_group_size', '') if arow else '',
            'overlap_strand_config': arow.get('overlap_strand_config', '') if arow else '',
            'overlap_n_pc':        arow.get('overlap_n_pc', '') if arow else '',
            'overlap_group_biotypes': arow.get('overlap_group_biotypes', '') if arow else '',
            'overlap_group_genes': arow.get('overlap_group_genes', '') if arow else '',
            'uf_overall':          uf_overall,
            'uf_3p250':            uf_3p250,
            'delta_3p250':         delta_3p250,
            'kallisto_tpm':        round(kal_tpm, 4) if kal_tpm is not None else '',
        }
        gene_rows.append(row)

    out_fields = ['enst_id', 'gene_symbol', 'transcript_len_nt', 'phenomenon',
                  'overlap_group_id', 'overlap_group_size', 'overlap_strand_config',
                  'overlap_n_pc', 'overlap_group_biotypes', 'overlap_group_genes',
                  'uf_overall', 'uf_3p250', 'delta_3p250', 'kallisto_tpm']
    with open(HERE / 'mane_overlap_genes.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields, delimiter='\t', lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(gene_rows)
    print(f'\nSaved mane_overlap_genes.tsv  ({len(gene_rows):,} MANE transcripts)')

    # ── Phenomenon UF summary ─────────────────────────────────────────────────
    phen_data = defaultdict(lambda: {'uf_overall': [], 'uf_3p250': [], 'delta_3p250': [],
                                      'kal_tpm': []})
    for r in gene_rows:
        p = r['phenomenon']
        if r['uf_overall'] != '':
            phen_data[p]['uf_overall'].append(float(r['uf_overall']))
        if r['uf_3p250'] != '':
            phen_data[p]['uf_3p250'].append(float(r['uf_3p250']))
        if r['delta_3p250'] != '':
            phen_data[p]['delta_3p250'].append(float(r['delta_3p250']))
        if r['kallisto_tpm'] != '':
            phen_data[p]['kal_tpm'].append(float(r['kallisto_tpm']))

    phen_rows = []
    for phen in ['unique', 'nested', 'genomic_overlap_antisense', 'genomic_overlap_sense',
                 'local_tandem_dup', 'readthrough', 'retrogene', 'pseudogene_shadow',
                 'distant_paralog', 'unknown']:
        d = phen_data.get(phen, {})
        ufo  = d.get('uf_overall', [])
        uf3  = d.get('uf_3p250', [])
        dlt  = d.get('delta_3p250', [])
        kal  = d.get('kal_tpm', [])
        if not ufo:
            continue
        phen_rows.append({
            'phenomenon':            phen,
            'n_genes':               len(ufo),
            'mean_uf_overall':       round(mean(ufo), 4),
            'median_uf_overall':     round(median(ufo), 4),
            'frac_uf_overall_lt0.5': round(sum(1 for u in ufo if u < 0.5) / len(ufo), 4),
            'mean_uf_3p250':         round(mean(uf3), 4) if uf3 else '',
            'median_uf_3p250':       round(median(uf3), 4) if uf3 else '',
            'frac_uf3p250_lt0.8':   round(sum(1 for u in uf3 if u < 0.8) / len(uf3), 4) if uf3 else '',
            'mean_delta_3p250':      round(mean(dlt), 4) if dlt else '',
            'mean_kallisto_tpm':     round(mean(kal), 2) if kal else '',
        })

    phen_fields = ['phenomenon', 'n_genes', 'mean_uf_overall', 'median_uf_overall',
                   'frac_uf_overall_lt0.5', 'mean_uf_3p250', 'median_uf_3p250',
                   'frac_uf3p250_lt0.8', 'mean_delta_3p250', 'mean_kallisto_tpm']
    with open(HERE / 'phenomenon_uf_summary.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=phen_fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(phen_rows)
    print(f'Saved phenomenon_uf_summary.tsv')

    # ── Overlap group summary ─────────────────────────────────────────────────
    # For genes that have an overlap_group_id, summarise read distribution
    # within each group (how is total kallisto TPM distributed across group members?)
    grp_genes = defaultdict(list)
    for r in gene_rows:
        ogid = r['overlap_group_id']
        if ogid and r['phenomenon'] not in ('unique', '', 'unknown'):
            grp_genes[ogid].append(r)

    grp_rows = []
    for ogid, members in sorted(grp_genes.items()):
        tpms  = [float(m['kallisto_tpm']) for m in members if m['kallisto_tpm'] != '']
        ufs   = [float(m['uf_overall'])   for m in members if m['uf_overall']   != '']
        total_tpm = sum(tpms)
        grp_rows.append({
            'overlap_group_id':      ogid,
            'n_members':             len(members),
            'gene_symbols':          ';'.join(m['gene_symbol'] for m in members),
            'phenomena':             ';'.join(sorted({m['phenomenon'] for m in members})),
            'strand_config':         members[0].get('overlap_strand_config', ''),
            'biotypes':              members[0].get('overlap_group_biotypes', ''),
            'n_with_uf':             len(ufs),
            'mean_uf_overall':       round(mean(ufs), 4) if ufs else '',
            'min_uf_overall':        round(min(ufs), 4) if ufs else '',
            'n_with_kallisto':       len(tpms),
            'total_kallisto_tpm':    round(total_tpm, 3),
            'max_tpm_fraction':      round(max(tpms) / total_tpm, 4) if total_tpm > 0 else '',
            'tpm_dominant_gene':     (
                members[tpms.index(max(tpms))]['gene_symbol']
                if tpms else ''
            ),
        })

    grp_fields = ['overlap_group_id', 'n_members', 'gene_symbols', 'phenomena',
                  'strand_config', 'biotypes', 'n_with_uf', 'mean_uf_overall',
                  'min_uf_overall', 'n_with_kallisto', 'total_kallisto_tpm',
                  'max_tpm_fraction', 'tpm_dominant_gene']
    with open(HERE / 'overlap_group_summary.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=grp_fields, delimiter='\t', lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(grp_rows)
    print(f'Saved overlap_group_summary.tsv  ({len(grp_rows)} overlap groups with MANE members)')

    # ── Console summary ───────────────────────────────────────────────────────
    print('\n═══ Summary ═══')
    print(f'{"Phenomenon":<30} {"n":>5}  {"mean UF_all":>11}  {"frac uf3p250<0.8":>16}')
    print('─' * 70)
    for r in phen_rows:
        print(f'{r["phenomenon"]:<30} {r["n_genes"]:>5}  '
              f'{r["mean_uf_overall"]:>11.4f}  '
              f'{str(r["frac_uf3p250_lt0.8"]):>16}')


if __name__ == '__main__':
    main()
