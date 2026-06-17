#!/usr/bin/env python3
"""
Gene-level summary of kallisto inflation analysis.

Reads the per-transcript inflation table (kallisto_inflation.tsv) and the
run summary (kallisto_run_summary.tsv) and produces:

  kallisto_genes_summary.tsv
    One row per transcript (ENST). For each gene: direction counts across all
    22 conditions, median TPM ratio (MANE/CDS), classification, and flags for
    families known to be problematic (RPL/RPS, HLA, HB*, HIST, etc.).

  kallisto_run_overview.tsv
    Compact one-row-per-condition table: read_source, mode, read_len,
    p_pseudoaligned (MANE), p_unique (MANE), p_pseudoaligned (CDS), p_unique (CDS),
    n_MANE_higher, n_CDS_higher, n_CDS_zero.

Output: results/kallisto_genes_summary.tsv
        results/kallisto_run_overview.tsv
"""

import csv
import statistics
from collections import defaultdict
from pathlib import Path

HERE    = Path(__file__).resolve().parent
RESULTS = HERE / 'results'

INFL_TSV = RESULTS / 'kallisto_inflation.tsv'
RUN_TSV  = RESULTS / 'kallisto_run_summary.tsv'

# Gene family prefixes to flag
FAMILIES = {
    'RPL':   'ribosomal_protein_L',
    'RPS':   'ribosomal_protein_S',
    'HLA-':  'HLA',
    'HBA':   'hemoglobin_alpha',
    'HBB':   'hemoglobin_beta',
    'HBD':   'hemoglobin_delta',
    'HIST':  'histone',
    'H2A':   'histone_H2A',
    'H2B':   'histone_H2B',
    'H3':    'histone_H3',
    'H4':    'histone_H4',
    'SNRP':  'snRNP',
    'UBB':   'ubiquitin',
    'UBC':   'ubiquitin',
    'CALM':  'calmodulin',
    'OR':    'olfactory_receptor',
    'TUBA':  'tubulin_alpha',
    'TUBB':  'tubulin_beta',
    'ACTB':  'actin',
    'ACTA':  'actin',
}


def gene_family(symbol):
    for prefix, name in FAMILIES.items():
        if symbol.startswith(prefix):
            return name
    return ''


# ── Load inflation table ────────────────────────────────────────────────────────

print("Loading inflation table ...", flush=True)
# {enst: {gene_symbol, [rows...]}}
gene_rows = defaultdict(list)
enst_to_gene = {}

with open(INFL_TSV) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        enst = row['transcript_id']
        gene_rows[enst].append(row)
        enst_to_gene[enst] = row['gene_symbol']

print(f"  {len(gene_rows):,} transcripts, "
      f"{len({v for v in enst_to_gene.values() if v}):,} unique gene symbols", flush=True)

# ── Gene-level summary ──────────────────────────────────────────────────────────

print("Computing gene-level summary ...", flush=True)

gene_out_rows = []

for enst, rows in sorted(gene_rows.items()):
    gsym = enst_to_gene.get(enst, '')

    n_total  = len(rows)
    dirs = [r['direction'] for r in rows]
    n_mane_higher = dirs.count('MANE_higher')
    n_cds_higher  = dirs.count('CDS_higher')
    n_equal       = dirs.count('equal')
    n_cds_zero    = dirs.count('CDS_zero')

    # Ratios only where both have signal
    ratios = []
    tpm_mane_vals = []
    tpm_cds_vals  = []
    for r in rows:
        try:
            ratio = float(r['ratio_mane_over_cds'])
            ratios.append(ratio)
        except (ValueError, TypeError):
            pass
        try:
            tpm_mane_vals.append(float(r['tpm_mane']))
            tpm_cds_vals.append(float(r['tpm_cds']))
        except (ValueError, TypeError):
            pass

    median_ratio = round(statistics.median(ratios), 4) if ratios else None
    median_tpm_mane = round(statistics.median(tpm_mane_vals), 4) if tpm_mane_vals else 0.0
    median_tpm_cds  = round(statistics.median(tpm_cds_vals),  4) if tpm_cds_vals  else 0.0

    # Classification
    frac_zero      = n_cds_zero / n_total
    frac_mane_hi   = n_mane_higher / n_total
    frac_cds_hi    = n_cds_higher / n_total

    if frac_zero >= 0.9:
        classification = 'no_CDS_signal'       # non-coding / absent from CDS ref
    elif frac_zero >= 0.5:
        classification = 'mostly_no_CDS_signal'
    elif frac_mane_hi >= 0.8 and median_ratio is not None and median_ratio > 2.0:
        classification = 'MANE_strongly_inflated'
    elif frac_mane_hi >= 0.6:
        classification = 'MANE_inflated'
    elif frac_cds_hi >= 0.8 and median_ratio is not None and median_ratio < 0.5:
        classification = 'CDS_strongly_inflated'
    elif frac_cds_hi >= 0.6:
        classification = 'CDS_inflated'
    elif frac_mane_hi >= 0.4 and frac_cds_hi >= 0.4:
        classification = 'mixed'
    else:
        classification = 'consistent'

    family = gene_family(gsym)

    gene_out_rows.append({
        'transcript_id':    enst,
        'gene_symbol':      gsym,
        'gene_family':      family,
        'n_conditions':     n_total,
        'n_MANE_higher':    n_mane_higher,
        'n_CDS_higher':     n_cds_higher,
        'n_equal':          n_equal,
        'n_CDS_zero':       n_cds_zero,
        'frac_MANE_higher': round(frac_mane_hi, 3),
        'frac_CDS_higher':  round(frac_cds_hi, 3),
        'frac_CDS_zero':    round(frac_zero, 3),
        'median_ratio_mane_over_cds': median_ratio if median_ratio is not None else 'NA',
        'median_tpm_mane':  median_tpm_mane,
        'median_tpm_cds':   median_tpm_cds,
        'classification':   classification,
    })

# Sort: no_CDS_signal first (non-coding), then by classification + ratio
cls_order = ['no_CDS_signal', 'mostly_no_CDS_signal', 'MANE_strongly_inflated',
             'MANE_inflated', 'CDS_strongly_inflated', 'CDS_inflated',
             'mixed', 'consistent']
cls_rank = {c: i for i, c in enumerate(cls_order)}
gene_out_rows.sort(key=lambda r: (
    cls_rank.get(r['classification'], 99),
    -(r['median_ratio_mane_over_cds'] if r['median_ratio_mane_over_cds'] != 'NA' else 0)
))

genes_path = RESULTS / 'kallisto_genes_summary.tsv'
gene_fields = ['transcript_id', 'gene_symbol', 'gene_family', 'n_conditions',
               'n_MANE_higher', 'n_CDS_higher', 'n_equal', 'n_CDS_zero',
               'frac_MANE_higher', 'frac_CDS_higher', 'frac_CDS_zero',
               'median_ratio_mane_over_cds', 'median_tpm_mane', 'median_tpm_cds',
               'classification']
with open(genes_path, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=gene_fields, delimiter='\t')
    w.writeheader()
    w.writerows(gene_out_rows)

print(f"  Saved → {genes_path}")

# ── Classification breakdown ───────────────────────────────────────────────────

from collections import Counter
cls_counts = Counter(r['classification'] for r in gene_out_rows)
fam_counts = Counter(r['gene_family'] for r in gene_out_rows if r['gene_family'])

print(f"\n{'═'*60}")
print(f"Classification summary  (n={len(gene_out_rows):,} transcripts)")
print(f"{'─'*60}")
for cls in cls_order:
    n = cls_counts.get(cls, 0)
    bar = '█' * (n // 200)
    print(f"  {cls:<30} {n:>6,}  {bar}")

print(f"\nFlagged gene families:")
for fam, n in sorted(fam_counts.items(), key=lambda x: -x[1])[:15]:
    # Count inflated in this family
    n_mane = sum(1 for r in gene_out_rows
                 if r['gene_family'] == fam
                 and 'MANE' in r['classification'])
    n_cds  = sum(1 for r in gene_out_rows
                 if r['gene_family'] == fam
                 and 'CDS' in r['classification']
                 and 'zero' not in r['classification'])
    print(f"  {fam:<30} n={n:>4}  MANE_inflated={n_mane}  CDS_inflated={n_cds}")

# ── Top MANE-inflated ───────────────────────────────────────────────────────────
print(f"\nTop 20 MANE-inflated genes (highest median ratio, CDS_zero excluded):")
print(f"  {'Gene':<12} {'Family':<25} {'median_ratio':>12}  {'class'}")
print(f"  {'─'*70}")
top_mane = [r for r in gene_out_rows
            if r['median_ratio_mane_over_cds'] != 'NA'
            and float(r['median_ratio_mane_over_cds']) > 1
            and r['frac_CDS_zero'] < 0.5]
top_mane.sort(key=lambda r: -float(r['median_ratio_mane_over_cds']))
for r in top_mane[:20]:
    print(f"  {r['gene_symbol']:<12} {r['gene_family']:<25} "
          f"{r['median_ratio_mane_over_cds']:>12}  {r['classification']}")

# ── Top CDS-inflated ────────────────────────────────────────────────────────────
print(f"\nTop 20 CDS-inflated genes (lowest median ratio, CDS_zero excluded):")
print(f"  {'Gene':<12} {'Family':<25} {'median_ratio':>12}  {'class'}")
print(f"  {'─'*70}")
top_cds = [r for r in gene_out_rows
           if r['median_ratio_mane_over_cds'] != 'NA'
           and float(r['median_ratio_mane_over_cds']) < 1
           and r['frac_CDS_zero'] < 0.5]
top_cds.sort(key=lambda r: float(r['median_ratio_mane_over_cds']))
for r in top_cds[:20]:
    print(f"  {r['gene_symbol']:<12} {r['gene_family']:<25} "
          f"{r['median_ratio_mane_over_cds']:>12}  {r['classification']}")

# ── Run overview table ──────────────────────────────────────────────────────────

print(f"\nBuilding run overview ...", flush=True)

run_rows = {}
with open(RUN_TSV) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        key = (row['read_source'], row['mode'], row['read_len_bp'])
        if key not in run_rows:
            run_rows[key] = {}
        idx = row['index']
        run_rows[key][idx] = row

# Count directions per condition from inflation table
dir_counts = defaultdict(Counter)
with open(INFL_TSV) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        key = (row['read_source'], row['mode'], row['read_len_bp'])
        dir_counts[key][row['direction']] += 1

overview_rows = []
ov_fields = ['read_source', 'mode', 'read_len_bp',
             'n_reads_processed',
             'MANE_p_pseudoaligned', 'MANE_p_unique',
             'CDS_p_pseudoaligned',  'CDS_p_unique',
             'n_MANE_higher', 'n_CDS_higher', 'n_equal', 'n_CDS_zero']

for key in sorted(run_rows):
    src, mode, rlen = key
    m = run_rows[key].get('MANE', {})
    c = run_rows[key].get('dedup_CDS', {})
    dc = dir_counts.get(key, {})
    overview_rows.append({
        'read_source':         src,
        'mode':                mode,
        'read_len_bp':         rlen,
        'n_reads_processed':   m.get('n_processed', ''),
        'MANE_p_pseudoaligned': m.get('p_pseudoaligned', ''),
        'MANE_p_unique':       m.get('p_unique', ''),
        'CDS_p_pseudoaligned': c.get('p_pseudoaligned', ''),
        'CDS_p_unique':        c.get('p_unique', ''),
        'n_MANE_higher':       dc.get('MANE_higher', 0),
        'n_CDS_higher':        dc.get('CDS_higher', 0),
        'n_equal':             dc.get('equal', 0),
        'n_CDS_zero':          dc.get('CDS_zero', 0),
    })

ov_path = RESULTS / 'kallisto_run_overview.tsv'
with open(ov_path, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=ov_fields, delimiter='\t')
    w.writeheader()
    w.writerows(overview_rows)

print(f"  Saved → {ov_path}")

# Print compact overview
print(f"\n{'─'*95}")
print(f"{'Read source':<15} {'Mo':>2} {'Len':>3}  "
      f"{'MANE %aln':>9} {'%uniq':>6}  "
      f"{'CDS %aln':>8} {'%uniq':>6}  "
      f"{'MANE>CDS':>9} {'CDS>MANE':>9} {'CDS=0':>7}")
print(f"{'─'*95}")
for r in overview_rows:
    print(f"{r['read_source']:<15} {r['mode']:>2} {r['read_len_bp']:>3}  "
          f"{r['MANE_p_pseudoaligned']:>9} {r['MANE_p_unique']:>6}  "
          f"{r['CDS_p_pseudoaligned']:>8} {r['CDS_p_unique']:>6}  "
          f"{r['n_MANE_higher']:>9,} {r['n_CDS_higher']:>9,} {r['n_CDS_zero']:>7,}")

print(f"\nSaved:\n  {genes_path}\n  {ov_path}")
