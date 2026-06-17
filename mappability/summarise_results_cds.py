#!/usr/bin/env python3
"""
Summarise v7 CDS-level uniqueness factor results.

v7 simulates reads from every position of each deduplicated CDS entry
(93,088 entries, NP_/XP_/YP_ accessions) and aligns them back to either
the full genome (8 conditions) or the CDS reference itself (2 self-map
conditions).  The unit is a CDS entry, not a transcript.

n_too_short = CDS entries absent from a given TSV because they are shorter
              than the simulated read length (no reads were generated).

Column breakdown per condition (all sum to n_total):
  n_too_short + f_eq_1 + f_0.10_to_1 + f_lt_0.10 = n_total
  where f_lt_0.10 includes f_eq_0 (reported separately)

Self-map note:
  In the selfmap conditions reads from a CDS entry align back to the same
  CDS reference. UF reflects within-CDS-reference redundancy: near-identical
  paralogs present as separate entries still cause multimapping. The selfmap
  UF is the ceiling — no aligner ambiguity from intronic sequences.

Output: results/cds_summary.tsv
"""

import csv
from pathlib import Path

HERE        = Path(__file__).resolve().parent
RESULTS_DIR = HERE / "results"

TSV_TO_SCRIPT = {
    'cds_uniqueness_factors_genome_SE.tsv':       'v7_a',
    'cds_uniqueness_factors_genome_SE_L75bp.tsv': 'v7_b',
    'cds_uniqueness_factors_genome_SE_L150bp.tsv':'v7_c',
    'cds_uniqueness_factors_genome_SE_L200bp.tsv':'v7_d',
    'cds_uniqueness_factors_genome_PE.tsv':        'v7_e',
    'cds_uniqueness_factors_genome_PE_L75bp.tsv':  'v7_f',
    'cds_uniqueness_factors_genome_PE_L150bp.tsv': 'v7_g',
    'cds_uniqueness_factors_genome_PE_L200bp.tsv': 'v7_h',
    'cds_uniqueness_factors_selfmap_SE.tsv':       'v7_i',
    'cds_uniqueness_factors_selfmap_PE.tsv':       'v7_j',
}

SCRIPT_ORDER = ['v7_a', 'v7_b', 'v7_c', 'v7_d',
                'v7_e', 'v7_f', 'v7_g', 'v7_h',
                'v7_i', 'v7_j']

# (reference, read_source, mode, read_len_bp)
METADATA = {
    'v7_a': ('full genome', 'CDS reads', 'SE', '100'),
    'v7_b': ('full genome', 'CDS reads', 'SE',  '75'),
    'v7_c': ('full genome', 'CDS reads', 'SE', '150'),
    'v7_d': ('full genome', 'CDS reads', 'SE', '200'),
    'v7_e': ('full genome', 'CDS reads', 'PE', '100'),
    'v7_f': ('full genome', 'CDS reads', 'PE',  '75'),
    'v7_g': ('full genome', 'CDS reads', 'PE', '150'),
    'v7_h': ('full genome', 'CDS reads', 'PE', '200'),
    'v7_i': ('CDS self',    'CDS reads', 'SE', '100'),
    'v7_j': ('CDS self',    'CDS reads', 'PE', '100'),
}

FIELDS = ['script', 'reference', 'read_source', 'mode', 'read_len_bp',
          'n_total', 'n_too_short', 'n_cds',
          'mean_uf', 'f_eq_1', 'f_ge_0.90', 'f_0.10_to_0.90', 'f_lt_0.10', 'f_eq_0',
          'total_positions', 'total_unique', 'total_multi', 'total_unmapped',
          'global_uf']


# ── Load all TSVs ──────────────────────────────────────────────────────────────

data    = {}   # {script: {cds_id: row_dict}}
all_ids = set()

for fname, script in TSV_TO_SCRIPT.items():
    d = {}
    path = RESULTS_DIR / fname
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            cid = row['cds_id']
            d[cid] = row
            all_ids.add(cid)
    data[script] = d
    print(f"  {script}  {fname}: {len(d):,} CDS entries")

n_total = len(all_ids)
print(f"\nTotal CDS entries (union across all TSVs): {n_total:,}\n")


# ── Compute per-condition stats ────────────────────────────────────────────────

out_path = RESULTS_DIR / 'cds_summary.tsv'
with open(out_path, 'w', newline='') as fh:
    writer = csv.writer(fh, delimiter='\t')
    writer.writerow(FIELDS)

    for script in SCRIPT_ORDER:
        d   = data[script]
        ufs = [float(row['uniqueness_factor']) for row in d.values()]
        n   = len(ufs)
        short = n_total - n

        f_eq_1  = sum(1 for f in ufs if f == 1.0)
        f_ge_90 = sum(1 for f in ufs if f >= 0.90)
        f_lt_10 = sum(1 for f in ufs if f <  0.10)
        f_eq_0  = sum(1 for f in ufs if f == 0.0)
        f_mid   = n - f_eq_1 - f_lt_10
        mean_uf = round(sum(ufs) / n, 4) if n else 0

        # Aggregate counts across all CDS entries for global UF
        tot_pos   = sum(int(row['n_positions'])    for row in d.values())
        tot_uniq  = sum(int(row['n_unique_back'])  for row in d.values())
        tot_multi = sum(int(row['n_multi_back'])   for row in d.values())
        tot_unmap = sum(int(row['n_unmapped'])     for row in d.values())
        global_uf = round(tot_uniq / tot_pos, 4) if tot_pos else 0

        ref, rsrc, mode, L = METADATA[script]
        writer.writerow([script, ref, rsrc, mode, L,
                         n_total, short, n,
                         mean_uf, f_eq_1, f_ge_90, f_mid, f_lt_10, f_eq_0,
                         tot_pos, tot_uniq, tot_multi, tot_unmap, global_uf])

        flag = ''
        check = short + f_eq_1 + f_mid + f_lt_10
        if check != n_total:
            flag = f'  *** CHECK FAILS ({check} != {n_total})'

        print(f"  {script} [{ref:11s} {mode} {L:>3}bp]  "
              f"too_short={short:4d}  n={n:,}  mean_uf={mean_uf:.4f}  "
              f"f=1:{f_eq_1:>6,}  f≥.90:{f_ge_90:>6,}  "
              f"f.10-.90:{f_mid:>6,}  f<.10:{f_lt_10:>6,}  f=0:{f_eq_0:>6,}  "
              f"global_uf={global_uf:.4f}{flag}")

print(f"\nSaved → {out_path}")
