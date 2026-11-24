#!/usr/bin/env python3
"""
06_unique_genes/analyze.py
Task: Catalogue truly unique genes — single copy, no overlap, unique mapping.

Context: The 'unique' phenomenon label in gene_overlap_annotation.tsv identifies
genes with no genomic overlap, no sequence cluster sharing, and no duplication
signal.  These are expected to have high UF (≈1.0) for both the overall transcript
and the 3' window.  They represent the reliable anchors for expression analysis —
well-behaved genes whose quantification is insensitive to aligner choice.

This script:
  1. Extracts protein-coding genes with phenomenon = 'unique'
  2. Verifies uniqueness through UF data (should be ≈1.0)
  3. Checks whether these genes are also conserved (GENCODE gene type, biotype)
  4. Compares STAR UF vs kallisto TPM — for truly unique genes both should agree
  5. Identifies any "unique" genes that unexpectedly show low UF (possible errors
     in the annotation or repeat elements inside the gene body)

Inputs (independent):
  ../sc3UTRs/results/gene_overlap_annotation.tsv
  ../sc3UTRs/results/sc3prime_v5g_MANE_RNA_PE_100.tsv
  ../kallisto_results/mappability_MANE_RNA_PE/kallisto_vs_MANE/abundance.tsv
  ../../Ref/MANE.GRCh38.v1.5.ensembl_rna.fna

Outputs:
  unique_genes.tsv             all 'unique' PC genes with UF and TPM
  unique_genes_anomalies.tsv   'unique' genes with unexpectedly low UF
  unique_vs_nonunique.tsv      summary comparison: unique vs each other phenomenon
"""

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean, median, stdev

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

LOW_UF_THRESHOLD = 0.8    # unexpected if 'unique' gene has UF < this


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
    # Both input result files use MANE RNA PE 100 bp (v5g for STAR UF,
    # mappability_MANE_RNA_PE for kallisto TPM).
    # The cross-phenomenon comparison in unique_vs_nonunique.tsv aggregates
    # all MANE genes from the same single run → no mixing, comparison valid.
    warn_if_mixing_runs([
        (SC3P_MANE, 100),   # STAR UF — MANE RNA PE 100 bp
        (KAL_ABUND, 100),   # kallisto TPM — MANE RNA PE 100 bp
    ])

    print('Loading gene overlap annotation …')
    annot = load_annotation(ANNOT)

    print('Loading MANE FASTA gene map …')
    fasta_map = load_mane_fasta_gene_map(MANE_FASTA)

    print('Loading STAR UF — MANE RNA PE 100 bp …')
    sc3p = load_sc3prime(SC3P_MANE)

    print('Loading kallisto TPM …')
    kal  = load_kallisto(KAL_ABUND)

    # Index annotation by gene_name for MANE lookup
    sym2annot = defaultdict(list)
    for gid, row in annot.items():
        sym2annot[row['gene_name']].append(row)

    # ── Build unique-gene rows ────────────────────────────────────────────────
    unique_rows = []
    anomaly_rows = []

    for enst, gsym in sorted(fasta_map.items(), key=lambda x: x[1]):
        arow_list = sym2annot.get(gsym, [])
        arow = None
        for c in arow_list:
            if c['phenomenon'] == 'unique' and c['gene_type'] == 'protein_coding':
                arow = c
                break
        if arow is None:
            continue   # not a 'unique' PC gene

        sc3row  = sc3p.get(enst, {})
        kal_tpm = kal.get(enst, None)

        ufo = sc3row.get('uf_overall', '')
        uf3 = sc3row.get('uf_3p250', '')
        dlt = sc3row.get('delta_3p250', '')

        row = {
            'enst_id':          enst,
            'gene_symbol':      gsym,
            'gene_id':          arow['gene_id'],
            'chrom':            arow.get('chrom', ''),
            'strand':           arow.get('strand', ''),
            'uf_overall':       ufo,
            'uf_3p250':         uf3,
            'delta_3p250':      dlt,
            'kallisto_tpm':     round(kal_tpm, 4) if kal_tpm is not None else '',
            'gencode_tags':     arow.get('gencode_tags', ''),
            'uf_status': (
                'high'    if ufo != '' and float(ufo) >= 0.9 else
                'medium'  if ufo != '' and float(ufo) >= LOW_UF_THRESHOLD else
                'low'     if ufo != '' and float(ufo) < LOW_UF_THRESHOLD else
                'missing'
            ),
        }
        unique_rows.append(row)

        if ufo != '' and float(ufo) < LOW_UF_THRESHOLD:
            anomaly_rows.append({**row, 'anomaly_reason':
                'repeat_in_body?' if float(ufo) < 0.5 else 'borderline'})

    u_fields = ['enst_id', 'gene_symbol', 'gene_id', 'chrom', 'strand',
                'uf_overall', 'uf_3p250', 'delta_3p250', 'kallisto_tpm',
                'gencode_tags', 'uf_status']
    with open(HERE / 'unique_genes.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=u_fields, delimiter='\t', lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(unique_rows)
    print(f'\nSaved unique_genes.tsv  ({len(unique_rows)} unique PC MANE genes)')

    a_fields = u_fields + ['anomaly_reason']
    with open(HERE / 'unique_genes_anomalies.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=a_fields, delimiter='\t', lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(anomaly_rows)
    print(f'Saved unique_genes_anomalies.tsv  ({len(anomaly_rows)} anomalous genes)')

    # ── Cross-phenomenon comparison ────────────────────────────────────────────
    # For all MANE genes (not just unique), group by phenomenon and compute UF stats
    phen_data = defaultdict(lambda: {'ufo': [], 'uf3': [], 'tpm': []})
    for enst, gsym in fasta_map.items():
        sc3row  = sc3p.get(enst, {})
        kal_tpm = kal.get(enst, None)
        ufo_s = sc3row.get('uf_overall', '')
        uf3_s = sc3row.get('uf_3p250', '')

        arow_list = sym2annot.get(gsym, [])
        phen = 'unknown'
        for c in arow_list:
            if c['gene_type'] == 'protein_coding':
                phen = c['phenomenon']
                break

        if ufo_s != '':
            phen_data[phen]['ufo'].append(float(ufo_s))
        if uf3_s != '':
            phen_data[phen]['uf3'].append(float(uf3_s))
        if kal_tpm is not None:
            phen_data[phen]['tpm'].append(kal_tpm)

    compare_rows = []
    for phen in ['unique', 'nested', 'genomic_overlap_antisense', 'genomic_overlap_sense',
                 'local_tandem_dup', 'readthrough', 'retrogene', 'pseudogene_shadow',
                 'distant_paralog']:
        d = phen_data.get(phen, {'ufo': [], 'uf3': [], 'tpm': []})
        ufo = d['ufo']
        uf3 = d['uf3']
        tpm = d['tpm']
        if not ufo:
            continue
        compare_rows.append({
            'phenomenon':         phen,
            'n_mane_genes':       len(ufo),
            'mean_uf_overall':    round(mean(ufo), 4),
            'median_uf_overall':  round(median(ufo), 4),
            'sd_uf_overall':      round(stdev(ufo), 4) if len(ufo) > 1 else 0,
            'frac_uf_ge_0.9':     round(sum(1 for u in ufo if u >= 0.9) / len(ufo), 4),
            'frac_uf_lt_0.8':     round(sum(1 for u in ufo if u < 0.8) / len(ufo), 4),
            'mean_uf_3p250':      round(mean(uf3), 4) if uf3 else '',
            'frac_uf3p250_lt0.8': round(sum(1 for u in uf3 if u < 0.8) / len(uf3), 4) if uf3 else '',
            'mean_kallisto_tpm':  round(mean(tpm), 2) if tpm else '',
        })

    comp_fields = ['phenomenon', 'n_mane_genes', 'mean_uf_overall', 'median_uf_overall',
                   'sd_uf_overall', 'frac_uf_ge_0.9', 'frac_uf_lt_0.8',
                   'mean_uf_3p250', 'frac_uf3p250_lt0.8', 'mean_kallisto_tpm']
    with open(HERE / 'unique_vs_nonunique.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=comp_fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(compare_rows)
    print(f'Saved unique_vs_nonunique.tsv')

    # ── Console ───────────────────────────────────────────────────────────────
    uf_high   = sum(1 for r in unique_rows if r['uf_status'] == 'high')
    uf_medium = sum(1 for r in unique_rows if r['uf_status'] == 'medium')
    uf_low    = sum(1 for r in unique_rows if r['uf_status'] == 'low')
    uf_miss   = sum(1 for r in unique_rows if r['uf_status'] == 'missing')
    n = len(unique_rows)
    print('\n═══ Summary ═══')
    print(f'Unique PC genes in MANE  : {n:,}')
    print(f'  UF ≥ 0.9  (high)       : {uf_high:,}  ({100*uf_high/n:.1f}%)')
    print(f'  UF 0.8-0.9 (medium)    : {uf_medium:,}  ({100*uf_medium/n:.1f}%)')
    print(f'  UF < 0.8  (anomalous)  : {uf_low:,}  ({100*uf_low/n:.1f}%)')
    print(f'  No UF data             : {uf_miss:,}')
    print(f'\nCross-phenomenon comparison:')
    print(f'{"Phenomenon":<30} {"n":>5}  {"mean UF":>7}  {"frac UF≥0.9":>11}  {"frac UF<0.8":>11}')
    print('─' * 70)
    for r in compare_rows:
        print(f'{r["phenomenon"]:<30} {r["n_mane_genes"]:>5}  '
              f'{r["mean_uf_overall"]:>7.4f}  '
              f'{r["frac_uf_ge_0.9"]:>11.4f}  '
              f'{r["frac_uf_lt_0.8"]:>11.4f}')


if __name__ == '__main__':
    main()
