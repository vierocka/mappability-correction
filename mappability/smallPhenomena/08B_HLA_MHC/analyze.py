#!/usr/bin/env python3
"""
08B_HLA_MHC/analyze.py
Case study: HLA/MHC antigen presentation pathway
Source: https://www.nature.com/articles/s41598-023-28167-1

Genes: HLA-A, HLA-B, HLA-C, HLA-E, HLA-F, HLA-G,
       HLA-DRA, HLA-DRB1, HLA-DQA1, HLA-DQB1, HLA-DQA2, HLA-DQB2,
       HLA-DPA1, HLA-DPB1, NLRC5, CIITA, PSME1, TAPBP, PSMB10, TAP2,
       HLA-DRB6, HLA-DMA, HLA-DOB, B2M, ERAP1, TAP1, PSME2, PSMB8,
       PSMB9, HLA-DRB5, HLA-DMB, HLA-DOA

The MHC locus (chr6p21.3) is one of the most gene-dense, polymorphic, and
repeat-rich regions in the human genome.  HLA class I and II genes have
thousands of alleles (allelic ambiguity), multiple pseudogene shadows,
and highly similar paralogous genes in the same locus (paralog ambiguity).

Some genes in this list (e.g., HLA-DRB6) are pseudogenes or non-MANE and
will not be found in the MANE RNA set — they are flagged in the output.

Read-length normalisation note (same as 08A)
----------------------------------------------
UF and uf_3pW are fractions [0,1] — cross-length comparison is directional.
TPM is per-million within each run.
est_counts must NOT be compared across runs.

Inputs (independent):
  ../../Ref/MANE.GRCh38.v1.5.ensembl_rna.fna
  ../sc3UTRs/results/gene_overlap_annotation.tsv
  ../sc3UTRs/results/sc3prime_*.tsv  (all MANE runs)
  ../kallisto_results/mappability_MANE_*/kallisto_vs_MANE/abundance.tsv

Outputs:
  gene_classification.tsv   one row per gene
  star_all_runs.tsv         one row per gene × STAR run
  kallisto_all_runs.tsv     one row per gene × kallisto run
  wide_uf_table.tsv         wide format for quick visual inspection
  mhc_locus_overview.tsv    genes sorted by chromosomal position, with overlap info
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
# Grouped by function for interpretation
TARGET_GENES = [
    # MHC class I (chr6p21.3) — allelic + paralog ambiguity
    'HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G',
    # MHC class II (chr6p21.3) — allelic + paralog ambiguity
    # HLA-DRB6 is a transcribed_unprocessed_pseudogene (GENCODE v49, no MANE transcript).
    # It is included because the Gong & Karchin paper listed it; it will appear as
    # "not in MANE" in gene_classification.tsv.  See HLA_DRB6_pseudogene_as_biomarker.txt
    # for the full analysis and connection to HLApseudo_MAFFT/classII_DRB alignment.
    'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DRB6',
    'HLA-DQA1', 'HLA-DQB1', 'HLA-DQA2', 'HLA-DQB2',
    'HLA-DPA1', 'HLA-DPB1',
    # HLA non-classical / DM / DO
    'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB',
    # Antigen presentation machinery (non-HLA)
    'B2M',     # beta-2-microglobulin, chr15 — unique gene
    'NLRC5',   # MHC-I transactivator, chr16 — NLRC/NLRP family
    'CIITA',   # MHC-II transactivator, chr16
    'TAPBP',   # Tapasin
    'TAP1',    'TAP2',   # Transporter associated with antigen processing
    'PSMB8',   'PSMB9',  'PSMB10',  # Immunoproteasome subunits
    'PSME1',   'PSME2',  # Proteasome activators
    'ERAP1',             # ER aminopeptidase
]

GENE_GROUPS = {
    'MHC_class_I':  ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G'],
    'MHC_class_II': ['HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DRB6',
                     'HLA-DQA1', 'HLA-DQB1', 'HLA-DQA2', 'HLA-DQB2',
                     'HLA-DPA1', 'HLA-DPB1', 'HLA-DMA', 'HLA-DMB',
                     'HLA-DOA', 'HLA-DOB'],
    'antigen_presentation': ['B2M', 'NLRC5', 'CIITA', 'TAPBP', 'TAP1', 'TAP2',
                             'PSMB8', 'PSMB9', 'PSMB10', 'PSME1', 'PSME2', 'ERAP1'],
}


def build_wide_table(target_genes, classifications, star_detail):
    mane_rna_labels = sorted(
        [lbl for lbl, meta in STAR_RUN_META.items()
         if meta[2] == 'RNA' and meta[3] == 'MANE'],
        key=lambda l: (STAR_RUN_META[l][1], STAR_RUN_META[l][0])
    )
    uf_by_gene_run = {}
    for r in star_detail:
        key = (r['gene_symbol'], r['run_label'])
        if key not in uf_by_gene_run:
            uf_by_gene_run[key] = r

    wide_rows = []
    for gsym in target_genes:
        cl = classifications.get(gsym, {})
        # Assign functional group
        group = next((g for g, genes in GENE_GROUPS.items() if gsym in genes), 'other')
        row = {
            'gene_symbol':    gsym,
            'functional_group': group,
            'ambiguity_types': cl.get('ambiguity_types', ''),
            'phenomenon':     cl.get('phenomenon', ''),
            'chrom':          cl.get('chrom', ''),
            'in_mane':        cl.get('in_mane', ''),
        }
        for lbl in mane_rna_labels:
            meta   = STAR_RUN_META[lbl]
            col    = f'uf_{meta[1]}_{meta[0]}bp'
            detail = uf_by_gene_run.get((gsym, lbl), {})
            row[col] = detail.get('uf_overall', '')
            row[f'uf3p250_{meta[1]}_{meta[0]}bp'] = detail.get('uf_3p250', '')
        wide_rows.append(row)

    wide_cols = (['gene_symbol', 'functional_group', 'ambiguity_types',
                  'phenomenon', 'chrom', 'in_mane'] +
                 [f'uf_{STAR_RUN_META[l][1]}_{STAR_RUN_META[l][0]}bp'
                  for l in mane_rna_labels] +
                 [f'uf3p250_{STAR_RUN_META[l][1]}_{STAR_RUN_META[l][0]}bp'
                  for l in mane_rna_labels])
    return wide_rows, wide_cols


def main():
    check_library_tsv(SC3UTR_DIR / 'sc3prime_v5g_MANE_RNA_PE_100.tsv', label_col='label')

    print(f'\nCase study 08B — HLA/MHC antigen presentation')
    print(f'{len(TARGET_GENES)} target genes\n')

    print('Building ENST → gene map from MANE FASTA …')
    enst_gene = build_enst_gene_map(MANE_FASTA)

    print('Loading gene overlap annotation …')
    gene_annot, _ = load_gene_annotation(ANNOT)

    print('Loading all MANE STAR runs …')
    star_runs = load_all_star_runs(SC3UTR_DIR, mane_only=True)
    print(f'  {len(star_runs)} runs: {sorted(star_runs.keys())}')

    print('Loading all MANE kallisto runs …')
    kal_runs = load_all_kallisto_runs(KAL_DIR)
    print(f'  {len(kal_runs)} runs')

    print('\nExtracting data for target genes …')
    classifications, star_detail, kal_detail = extract_gene_rows(
        TARGET_GENES, enst_gene, gene_annot, star_runs, kal_runs
    )

    gene_order = {g: i for i, g in enumerate(TARGET_GENES)}

    # ── gene_classification.tsv ───────────────────────────────────────────────
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

    # ── star_all_runs.tsv ─────────────────────────────────────────────────────
    star_fields = ['gene_symbol', 'enst_id', 'run_label', 'read_len_bp', 'mode',
                   'read_source', 'reference', 'transcript_len_nt',
                   'uf_overall', 'uf_3p250', 'delta_3p250',
                   'uf_3p400', 'delta_3p400', 'uf_3p500', 'delta_3p500',
                   'n_positions', 'n_unique_3p250', 'ambiguity_types']
    with open(HERE / 'star_all_runs.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=star_fields, delimiter='\t',
                           lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        for r in sorted(star_detail,
                        key=lambda x: (gene_order.get(x['gene_symbol'], 99),
                                       x['run_label'])):
            w.writerow(r)
    print(f'Saved star_all_runs.tsv  ({len(star_detail)} rows)')

    # ── kallisto_all_runs.tsv ─────────────────────────────────────────────────
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

    # ── wide_uf_table.tsv ─────────────────────────────────────────────────────
    wide_rows, wide_cols = build_wide_table(TARGET_GENES, classifications, star_detail)
    with open(HERE / 'wide_uf_table.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=wide_cols, delimiter='\t',
                           lineterminator='\n', extrasaction='ignore', restval='')
        w.writeheader()
        w.writerows(wide_rows)
    print(f'Saved wide_uf_table.tsv')

    # ── mhc_locus_overview.tsv — chr6 genes sorted by position ───────────────
    locus_rows = []
    for gsym in TARGET_GENES:
        cl = classifications.get(gsym, {})
        if cl.get('chrom') != 'chr6':
            continue
        # Get position from annotation
        arow = gene_annot.get(gsym, {})
        locus_rows.append({
            'gene_symbol':    gsym,
            'chrom':          cl.get('chrom', ''),
            'start':          arow.get('start', ''),
            'end':            arow.get('end', ''),
            'strand':         cl.get('strand', ''),
            'gene_type':      cl.get('gene_type', ''),
            'phenomenon':     cl.get('phenomenon', ''),
            'ambiguity_types': cl.get('ambiguity_types', ''),
            'overlap_group_size': cl.get('overlap_group_size', ''),
            'overlap_group_genes': cl.get('overlap_group_genes', ''),
        })
    # Sort by start position
    locus_rows.sort(key=lambda r: int(r['start']) if r['start'] else 0)
    locus_fields = ['gene_symbol', 'chrom', 'start', 'end', 'strand',
                    'gene_type', 'phenomenon', 'ambiguity_types',
                    'overlap_group_size', 'overlap_group_genes']
    with open(HERE / 'mhc_locus_overview.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=locus_fields, delimiter='\t',
                           lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        w.writerows(locus_rows)
    print(f'Saved mhc_locus_overview.tsv  ({len(locus_rows)} chr6 genes)')

    # ── Console summary ───────────────────────────────────────────────────────
    print('\n═══ Gene classification ═══')
    print(f'{"Gene":<12} {"Group":<22} {"Chrom":<6} {"in_MANE":<8} '
          f'{"Phenomenon":<28} {"Ambiguity types"}')
    print('─' * 110)
    for gsym in TARGET_GENES:
        cl    = classifications.get(gsym, {})
        group = next((g for g, genes in GENE_GROUPS.items() if gsym in genes), 'other')
        print(f'{gsym:<12} {group:<22} {cl.get("chrom","?"):<6} '
              f'{str(cl.get("in_mane","?")):<8} '
              f'{cl.get("phenomenon","?"):<28} {cl.get("ambiguity_types","?")}')

    # Warn about genes not in MANE
    missing = [g for g, cl in classifications.items() if not cl.get('in_mane')]
    if missing:
        print(f'\nGenes not in MANE (no UF data): {", ".join(missing)}')
        print('  These may be pseudogenes, non-canonical, or absent from MANE v1.5.')


if __name__ == '__main__':
    main()
