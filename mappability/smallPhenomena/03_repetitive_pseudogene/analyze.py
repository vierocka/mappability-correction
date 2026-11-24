#!/usr/bin/env python3
"""
03_repetitive_pseudogene/analyze.py
Task: Genes in highly repetitive regions with many pseudogenes.

Context: Many protein-coding genes have one or more processed or unprocessed
pseudogene copies elsewhere in the genome.  Reads from the coding gene can map
to the pseudogene locus (and vice versa), inflating or deflating quantification.
The CDS is usually the most conserved part; UTRs can help to distinguish the
functional copy from its pseudogene shadow.

This script:
  1. Loads GENCODE overlap annotation; identifies genes in 'pseudogene_shadow'
     phenomenon, plus protein-coding genes in groups that contain pseudogenes
     (overlap_group_biotypes contains 'pseudogene').
  2. For MANE transcripts in those groups, extracts:
       - uf_overall    (from MANE RNA PE 100 bp run)
       - uf_3p250      (3' end uniqueness — UTR region)
       - delta_3p250   (improvement from UTR over CDS+body average)
  3. Loads the CDS self-map uniqueness (STAR: CDS reads → dedup_cds) to show
     how non-unique the CDS portion is.
  4. Contrasts CDS uniqueness vs 3' UTR uniqueness: UTR benefit.

Inputs (independent):
  ../sc3UTRs/results/gene_overlap_annotation.tsv
  ../sc3UTRs/results/sc3prime_v5g_MANE_RNA_PE_100.tsv
  ../../Ref/MANE.GRCh38.v1.5.ensembl_rna.fna
  ../results/id_conversion_table.tsv
  ../sc3UTRs/results/sc3prime_v2e_dedupcds_RNA_PE_100.tsv   (CDS-only MANE run)

Outputs:
  pseudogene_genes.tsv       protein-coding genes that co-occur with pseudogenes
  pseudogene_utr_benefit.tsv per-gene: CDS-UF vs 3'UTR-UF delta (UTR rescue?)
  pseudogene_summary.tsv     group-level stats
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

ANNOT     = MAPDIR / 'sc3UTRs' / 'results' / 'gene_overlap_annotation.tsv'
SC3P_MANE = MAPDIR / 'sc3UTRs' / 'results' / 'sc3prime_v5g_MANE_RNA_PE_100.tsv'
SC3P_CDS  = MAPDIR / 'sc3UTRs' / 'results' / 'sc3prime_v2e_dedupcds_RNA_PE_100.tsv'
ID_TABLE  = MAPDIR / 'results' / 'id_conversion_table.tsv'
MANE_FASTA = REFDIR / 'MANE.GRCh38.v1.5.ensembl_rna.fna'

PSEUDOGENE_BIOTYPES = {
    'processed_pseudogene', 'unprocessed_pseudogene',
    'transcribed_unprocessed_pseudogene', 'transcribed_processed_pseudogene',
    'rRNA_pseudogene', 'transcribed_unitary_pseudogene', 'unitary_pseudogene',
    'IG_V_pseudogene', 'IG_C_pseudogene', 'IG_D_pseudogene', 'IG_J_pseudogene',
    'TR_V_pseudogene', 'TR_J_pseudogene', 'polymorphic_pseudogene',
    'pseudogene',
}


def has_pseudogene_biotype(biotypes_str):
    return any(b in PSEUDOGENE_BIOTYPES for b in biotypes_str.split(';'))


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


def load_id_table(path):
    enst2gene = {}
    gene2enst = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            enst = row.get('ensembl_transcript_id', '').split('.')[0]
            gsym = row.get('gene_symbol', '')
            if enst and gsym:
                enst2gene[enst] = gsym
                gene2enst[gsym] = enst
    return enst2gene, gene2enst


def load_mane_fasta_gene_map(path):
    """Returns {ENST_no_version: gene_symbol}."""
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
    # v5g (MANE RNA PE 100 bp) and v2e (dedup CDS RNA PE 100 bp) both use 100 bp PE.
    # Comparing UF across these two runs is valid (same read length).
    # However, note that v5g maps to the full MANE index (UTRs included) while
    # v2e maps to the dedup CDS index (CDS only) — so UF values reflect
    # different reference spaces, not just different library sizes.
    warn_if_mixing_runs([
        (SC3P_MANE, 100),   # v5g MANE RNA PE 100 bp (full transcriptome reference)
        (SC3P_CDS,  100),   # v2e dedup CDS RNA PE 100 bp (CDS-only reference)
    ])

    print('Loading gene overlap annotation …')
    annot = load_annotation(ANNOT)

    print('Loading ID table …')
    enst2gene, gene2enst = load_id_table(ID_TABLE)

    print('Loading MANE FASTA gene map …')
    fasta_enst2gene = load_mane_fasta_gene_map(MANE_FASTA)

    print('Loading STAR UF — MANE RNA PE 100 bp …')
    sc3p_mane = load_sc3prime(SC3P_MANE)

    print('Loading STAR UF — dedup CDS RNA PE 100 bp …')
    sc3p_cds = load_sc3prime(SC3P_CDS)

    # ── Identify protein-coding MANE genes in pseudogene-containing groups ────
    # Criterion 1: phenomenon == 'pseudogene_shadow'
    # Criterion 2: protein-coding gene in a group where biotypes contain pseudogene
    pseudo_groups = set()
    pc_in_pseudo  = {}   # ensg → annot row (protein-coding genes near pseudogenes)
    for gid, row in annot.items():
        biotypes = row.get('overlap_group_biotypes', '')
        if row['gene_type'] == 'protein_coding' and has_pseudogene_biotype(biotypes):
            pseudo_groups.add(row.get('overlap_group_id', ''))
            pc_in_pseudo[gid] = row
        elif row['phenomenon'] == 'pseudogene_shadow':
            pc_in_pseudo[gid] = row

    print(f'  {len(pc_in_pseudo):,} protein-coding genes co-located with pseudogenes')
    print(f'  {len(pseudo_groups):,} unique overlap groups')

    # ── Build output rows ─────────────────────────────────────────────────────
    out_rows = []
    utr_rows = []

    for enst, gsym in sorted(fasta_enst2gene.items(), key=lambda x: x[1]):
        # Find annotation by gene symbol
        arow = None
        for gid, row in annot.items():
            if row['gene_name'] == gsym and row['gene_type'] == 'protein_coding':
                arow = row
                break
        if arow is None:
            continue
        gid = arow['gene_id']
        if gid not in pc_in_pseudo:
            continue

        mane_row = sc3p_mane.get(enst, {})
        cds_row  = sc3p_cds.get(enst, {})

        uf_mane    = mane_row.get('uf_overall', '')
        uf_3p250   = mane_row.get('uf_3p250', '')
        delta_3p   = mane_row.get('delta_3p250', '')
        uf_cds_run = cds_row.get('uf_overall', '')   # same transcript but CDS index run

        row_out = {
            'enst_id':              enst,
            'gene_symbol':          gsym,
            'gene_id':              gid,
            'phenomenon':           arow['phenomenon'],
            'overlap_group_id':     arow.get('overlap_group_id', ''),
            'overlap_group_size':   arow.get('overlap_group_size', ''),
            'overlap_group_biotypes': arow.get('overlap_group_biotypes', ''),
            'n_pseudogenes_in_group': arow.get('overlap_n_pseudogene', ''),
            'uf_overall_mane_run':  uf_mane,
            'uf_3p250_mane_run':    uf_3p250,
            'delta_3p250':          delta_3p,
            'uf_overall_cds_run':   uf_cds_run,
            'utr_rescues':          (
                'yes' if uf_3p250 != '' and uf_mane != '' and
                float(uf_3p250) > float(uf_mane) + 0.1 else
                'partial' if uf_3p250 != '' and uf_mane != '' and
                float(uf_3p250) > float(uf_mane) else
                'no'
            ) if uf_3p250 != '' and uf_mane != '' else '',
        }
        out_rows.append(row_out)

        # UTR benefit row
        if uf_mane != '' and uf_3p250 != '':
            utr_rows.append({
                'gene_symbol':   gsym,
                'uf_cds_body':   uf_mane,     # overall UF dominated by CDS exons
                'uf_3p250':      uf_3p250,    # uniqueness of 3' end (UTR region)
                'delta_3p250':   delta_3p,    # positive = UTR helps, negative = UTR hurts
                'utr_benefit':   round(float(uf_3p250) - float(uf_mane), 4),
                'n_pseudo_in_group': arow.get('overlap_n_pseudogene', ''),
            })

    out_fields = ['enst_id', 'gene_symbol', 'gene_id', 'phenomenon',
                  'overlap_group_id', 'overlap_group_size', 'overlap_group_biotypes',
                  'n_pseudogenes_in_group', 'uf_overall_mane_run', 'uf_3p250_mane_run',
                  'delta_3p250', 'uf_overall_cds_run', 'utr_rescues']
    with open(HERE / 'pseudogene_genes.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields, delimiter='\t', lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(out_rows)
    print(f'\nSaved pseudogene_genes.tsv  ({len(out_rows)} genes)')

    utr_fields = ['gene_symbol', 'uf_cds_body', 'uf_3p250', 'delta_3p250',
                  'utr_benefit', 'n_pseudo_in_group']
    with open(HERE / 'pseudogene_utr_benefit.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=utr_fields, delimiter='\t', lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(utr_rows)
    print(f'Saved pseudogene_utr_benefit.tsv  ({len(utr_rows)} genes with UF data)')

    # ── Summary ───────────────────────────────────────────────────────────────
    utr_benefits = [float(r['utr_benefit']) for r in utr_rows]
    rescued = sum(1 for r in out_rows if r.get('utr_rescues') == 'yes')
    partial = sum(1 for r in out_rows if r.get('utr_rescues') == 'partial')
    no_resc = sum(1 for r in out_rows if r.get('utr_rescues') == 'no')

    phen_counts = defaultdict(int)
    for r in out_rows:
        phen_counts[r['phenomenon']] += 1

    summary_rows = [
        {'metric': 'n_genes_near_pseudogenes', 'value': len(out_rows)},
        {'metric': 'n_pseudogene_shadow',       'value': phen_counts.get('pseudogene_shadow', 0)},
        {'metric': 'n_other_near_pseudo',        'value': len(out_rows) - phen_counts.get('pseudogene_shadow', 0)},
        {'metric': 'mean_utr_benefit',           'value': round(mean(utr_benefits), 4) if utr_benefits else 'NA'},
        {'metric': 'frac_utr_rescues_yes',       'value': round(rescued / len(out_rows), 4) if out_rows else 'NA'},
        {'metric': 'frac_utr_rescues_partial',   'value': round(partial / len(out_rows), 4) if out_rows else 'NA'},
        {'metric': 'frac_utr_no_rescue',         'value': round(no_resc / len(out_rows), 4) if out_rows else 'NA'},
    ]
    with open(HERE / 'pseudogene_summary.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['metric', 'value'], delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(summary_rows)
    print(f'Saved pseudogene_summary.tsv')

    print('\n═══ Summary ═══')
    for r in summary_rows:
        print(f'  {r["metric"]:<35}: {r["value"]}')


if __name__ == '__main__':
    main()
