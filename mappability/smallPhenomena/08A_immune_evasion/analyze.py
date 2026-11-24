#!/usr/bin/env python3
"""
08A_immune_evasion/analyze.py
Case study: Immune evasion paper gene panel
Source: https://www.nature.com/articles/s41598-022-15482-2

Genes: GBP1, HLA-B, ELF4, HLA-E, IFITM2, TNFRSF14, CD274, AIM2, CFB, RHOG, HORMAD1

For each gene:
  - Classify by ambiguity type (genomic, isoform, paralog, pseudogene, allelic)
  - Pull STAR UF across ALL MANE runs (75/100/150/200 bp, PE/SE, RNA/exflank, proper-pair)
  - Pull kallisto TPM across ALL MANE runs (same read lengths and modes)
  - Compare behaviour across read lengths (key: longer reads → higher UF for ambiguous genes)

Read-length normalisation note
-------------------------------
UF (uniqueness factor) and uf_3pW are fractions [0,1] — safe to compare
directly across read lengths as RELATIVE measures.
TPM is per-million within each run — safe for within-run gene comparisons,
and the RATIO of a gene's TPM to its expected value can be compared across
runs as a fraction.
est_counts are NOT normalised by library size — do not compare raw counts
across runs (column is labelled accordingly in output).

Inputs (independent — no dependency on other smallPhenomena tasks):
  ../../Ref/MANE.GRCh38.v1.5.ensembl_rna.fna
  ../sc3UTRs/results/gene_overlap_annotation.tsv
  ../sc3UTRs/results/sc3prime_*.tsv   (all MANE runs)
  ../kallisto_results/mappability_MANE_*/kallisto_vs_MANE/abundance.tsv

Outputs:
  gene_classification.tsv   one row per gene: ambiguity types + annotation
  star_all_runs.tsv         one row per gene × STAR run: all UF metrics
  kallisto_all_runs.tsv     one row per gene × kallisto run: TPM
  wide_uf_table.tsv         wide format: gene rows, one col per run (for quick scan)
"""

import csv
import sys
from pathlib import Path

HERE   = Path(__file__).resolve().parent
ROOT   = HERE.parent.parent.parent
REFDIR = ROOT / 'Ref'
MAPDIR = ROOT / 'mappability'

sys.path.insert(0, str(HERE.parent))
from lib_check import check_library_tsv
from case_study_lib import (
    build_enst_gene_map, load_gene_annotation,
    load_all_star_runs, load_all_kallisto_runs,
    extract_gene_rows, STAR_RUN_META, KAL_RUN_META,
)

MANE_FASTA = REFDIR / 'MANE.GRCh38.v1.5.ensembl_rna.fna'
ANNOT      = MAPDIR / 'sc3UTRs' / 'results' / 'gene_overlap_annotation.tsv'
SC3UTR_DIR = MAPDIR / 'sc3UTRs' / 'results'
KAL_DIR    = MAPDIR / 'kallisto_results'

# ── Gene list from the paper ───────────────────────────────────────────────────
TARGET_GENES = [
    'GBP1',      # Guanylate binding protein 1 — GBP family paralog
    'HLA-B',     # MHC class I — allelic + paralog ambiguity
    'ELF4',      # ETS transcription factor
    'HLA-E',     # MHC class I (non-classical) — allelic + paralog
    'IFITM2',    # Interferon transmembrane protein — IFITM family
    'TNFRSF14',  # TNF receptor superfamily — large paralog family
    'CD274',     # PD-L1 — immune checkpoint
    'AIM2',      # ALR/HIN domain — inflammasome
    'CFB',       # Complement factor B — near MHC locus
    'RHOG',      # Rho GTPase — large paralog family
    'HORMAD1',   # HORMA domain protein
]


def build_wide_table(classifications, star_detail):
    """
    Pivot star_detail into a wide table:
    rows = genes, columns = run_label (uf_overall) for MANE RNA runs only.
    Helps quick visual comparison across read lengths.
    """
    # Ordered MANE RNA runs by (mode, read_len)
    mane_rna_labels = sorted(
        [lbl for lbl, meta in STAR_RUN_META.items()
         if meta[2] == 'RNA' and meta[3] == 'MANE'],
        key=lambda l: (STAR_RUN_META[l][1], STAR_RUN_META[l][0])
    )

    # Index star_detail by (gene, label)
    uf_by_gene_run = {}
    for r in star_detail:
        key = (r['gene_symbol'], r['run_label'])
        # If multiple ENST per gene (shouldn't happen for MANE), take first
        if key not in uf_by_gene_run:
            uf_by_gene_run[key] = r

    wide_rows = []
    for gsym in TARGET_GENES:
        cl = classifications.get(gsym, {})
        row = {
            'gene_symbol':    gsym,
            'ambiguity_types': cl.get('ambiguity_types', ''),
            'phenomenon':     cl.get('phenomenon', ''),
            'chrom':          cl.get('chrom', ''),
        }
        for lbl in mane_rna_labels:
            meta   = STAR_RUN_META[lbl]
            col    = f'uf_{meta[1]}_{meta[0]}bp'   # e.g. uf_PE_100bp
            detail = uf_by_gene_run.get((gsym, lbl), {})
            row[col] = detail.get('uf_overall', '')
            # Also 3' window
            row[f'uf3p250_{meta[1]}_{meta[0]}bp'] = detail.get('uf_3p250', '')
        wide_rows.append(row)

    wide_cols = (['gene_symbol', 'ambiguity_types', 'phenomenon', 'chrom'] +
                 [f'uf_{STAR_RUN_META[l][1]}_{STAR_RUN_META[l][0]}bp'
                  for l in mane_rna_labels] +
                 [f'uf3p250_{STAR_RUN_META[l][1]}_{STAR_RUN_META[l][0]}bp'
                  for l in mane_rna_labels])
    return wide_rows, wide_cols


def main():
    # ── Library consistency check ─────────────────────────────────────────────
    # We load ALL sc3prime files (multiple read lengths and modes).
    # UF/uf_3pW are fractions — cross-length comparison is directional (longer
    # reads generally give higher UF), not a normalisation failure.
    # TPM values from kallisto are per-million within each run.
    # est_counts are explicitly labelled as non-comparable across runs.
    check_library_tsv(SC3UTR_DIR / 'sc3prime_v5g_MANE_RNA_PE_100.tsv', label_col='label')

    print(f'\nCase study 08A — Immune evasion paper')
    print(f'Genes: {", ".join(TARGET_GENES)}\n')

    print('Building ENST → gene map from MANE FASTA …')
    enst_gene = build_enst_gene_map(MANE_FASTA)
    print(f'  {len(enst_gene):,} MANE transcripts')

    print('Loading gene overlap annotation …')
    gene_annot, _ = load_gene_annotation(ANNOT)

    print('Loading all MANE STAR runs …')
    star_runs = load_all_star_runs(SC3UTR_DIR, mane_only=True)
    print(f'  {len(star_runs)} MANE STAR runs loaded: {sorted(star_runs.keys())}')

    print('Loading all MANE kallisto runs …')
    kal_runs = load_all_kallisto_runs(KAL_DIR)
    print(f'  {len(kal_runs)} MANE kallisto runs loaded')

    print('\nExtracting data for target genes …')
    classifications, star_detail, kal_detail = extract_gene_rows(
        TARGET_GENES, enst_gene, gene_annot, star_runs, kal_runs
    )

    # ── Save classification table ─────────────────────────────────────────────
    cl_fields = ['gene_symbol', 'in_mane', 'n_enst_ids', 'enst_ids',
                 'chrom', 'strand', 'gene_type', 'phenomenon',
                 'ambiguity_types', 'is_HLA',
                 'overlap_group_size', 'overlap_group_biotypes', 'overlap_group_genes',
                 'seq_cluster_n_pc', 'seq_cluster_max_span_mb']
    with open(HERE / 'gene_classification.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cl_fields, delimiter='\t',
                           lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        for gsym in TARGET_GENES:
            if gsym in classifications:
                w.writerow(classifications[gsym])
    print(f'Saved gene_classification.tsv')

    # ── Save STAR all-runs table ──────────────────────────────────────────────
    star_fields = ['gene_symbol', 'enst_id', 'run_label', 'read_len_bp', 'mode',
                   'read_source', 'reference', 'transcript_len_nt',
                   'uf_overall', 'uf_3p250', 'delta_3p250',
                   'uf_3p400', 'delta_3p400', 'uf_3p500', 'delta_3p500',
                   'n_positions', 'n_unique_3p250', 'ambiguity_types']
    with open(HERE / 'star_all_runs.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=star_fields, delimiter='\t',
                           lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        # Sort by gene (in TARGET_GENES order), then run_label
        gene_order = {g: i for i, g in enumerate(TARGET_GENES)}
        for r in sorted(star_detail,
                        key=lambda x: (gene_order.get(x['gene_symbol'], 99),
                                       x['run_label'])):
            w.writerow(r)
    print(f'Saved star_all_runs.tsv  ({len(star_detail)} rows)')

    # ── Save kallisto all-runs table ──────────────────────────────────────────
    kal_fields = ['gene_symbol', 'enst_id', 'run_name', 'read_len_bp', 'mode',
                  'read_source', 'reference', 'tpm',
                  'est_counts_DO_NOT_COMPARE_ACROSS_RUNS', 'ambiguity_types']
    with open(HERE / 'kallisto_all_runs.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=kal_fields, delimiter='\t',
                           lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        for r in sorted(kal_detail,
                        key=lambda x: (gene_order.get(x['gene_symbol'], 99),
                                       x['run_name'])):
            w.writerow(r)
    print(f'Saved kallisto_all_runs.tsv  ({len(kal_detail)} rows)')

    # ── Save wide UF table ────────────────────────────────────────────────────
    wide_rows, wide_cols = build_wide_table(classifications, star_detail)
    with open(HERE / 'wide_uf_table.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=wide_cols, delimiter='\t',
                           lineterminator='\n', extrasaction='ignore', restval='')
        w.writeheader()
        w.writerows(wide_rows)
    print(f'Saved wide_uf_table.tsv  ({len(wide_cols)-4} UF columns)')

    # ── Console classification summary ────────────────────────────────────────
    print('\n═══ Gene classification ═══')
    print(f'{"Gene":<12} {"Chrom":<6} {"Phenomenon":<30} {"Ambiguity types"}')
    print('─' * 90)
    for gsym in TARGET_GENES:
        cl = classifications.get(gsym, {})
        print(f'{gsym:<12} {cl.get("chrom","?"):<6} '
              f'{cl.get("phenomenon","?"):<30} {cl.get("ambiguity_types","?")}')


if __name__ == '__main__':
    main()
