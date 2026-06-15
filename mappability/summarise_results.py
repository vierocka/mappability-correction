#!/usr/bin/env python3
"""
Recompute results_summary.tsv from the per-transcript uniqueness TSVs.

n_total  = union of all transcript IDs across all 13 TSVs (19,437)
tr_too_short = transcripts absent from a given TSV because they are shorter
               than the simulated read length (no reads could be generated)

Column breakdown per script (all sum to n_total):
  tr_too_short + f_eq_1 + f_0.10_to_1 + f_lt_0.10 = n_total
  where f_lt_0.10 includes f_eq_0 (reported separately)

Output: results/results_summary.tsv
"""

import csv
from pathlib import Path

HERE        = Path(__file__).resolve().parent
RESULTS_DIR = HERE / "results"

TSV_TO_SCRIPT = {
    'transcript_uniqueness_factors_genomic_L100bp.tsv':                  'v1_a',
    'transcript_uniqueness_factors_genomic_SE_L100bp.tsv':               'v1_b',
    'transcript_uniqueness_factors_genomic_RNA_SE_L100bp.tsv':           'v1_c',
    'transcript_uniqueness_factors_genomic_L75bp.tsv':                   'v3_a',
    'transcript_uniqueness_factors_genomic_L150bp.tsv':                  'v4_a',
    'transcript_uniqueness_factors_dedup_cds_PE_flank_L100bp.tsv':       'v2_a',
    'transcript_uniqueness_factors_dedup_cds_SE_flank_L100bp.tsv':       'v2_b',
    'transcript_uniqueness_factors_dedup_cds_L100bp.tsv':                'v2_c',
    'transcript_uniqueness_factors_dedup_cds_alignIntronMax1_L100bp.tsv':'v2_d',
    'transcript_uniqueness_factors_MANE_RNA_SE_L100bp.tsv':              'v5_a',
    'transcript_uniqueness_factors_MANE_exonflank_PE_L100bp.tsv':        'v5_b',
    'transcript_uniqueness_factors_MANE_RNA_SE_L75bp.tsv':               'v5_c',
    'transcript_uniqueness_factors_MANE_RNA_SE_L200bp.tsv':              'v5_d',
}

SCRIPT_ORDER = ['v1_a','v1_b','v1_c','v3_a','v4_a',
                'v2_a','v2_b','v2_c','v2_d',
                'v5_a','v5_b','v5_c','v5_d']

METADATA = {
    'v1_a': ('full genome', 'exon-flank', 'PE',  '100', 'no'),
    'v1_b': ('full genome', 'exon-flank', 'SE',  '100', 'no'),
    'v1_c': ('full genome', 'MANE RNA',   'SE',  '100', 'no'),
    'v3_a': ('full genome', 'exon-flank', 'PE',   '75', 'no'),
    'v4_a': ('full genome', 'exon-flank', 'PE',  '150', 'no'),
    'v2_a': ('dedup CDS',   'exon-flank', 'PE',  '100', 'no'),
    'v2_b': ('dedup CDS',   'exon-flank', 'SE',  '100', 'no'),
    'v2_c': ('dedup CDS',   'MANE RNA',   'SE',  '100', 'no'),
    'v2_d': ('dedup CDS',   'MANE RNA',   'SE',  '100', 'yes'),
    'v5_a': ('MANE ref',    'MANE RNA',   'SE',  '100', 'no'),
    'v5_b': ('MANE ref',    'exon-flank', 'PE',  '100', 'no'),
    'v5_c': ('MANE ref',    'MANE RNA',   'SE',   '75', 'no'),
    'v5_d': ('MANE ref',    'MANE RNA',   'SE',  '200', 'no'),
}

FIELDS = ['script', 'reference', 'read_source', 'mode', 'read_len_bp', 'alignIntronMax1',
          'n_total', 'tr_too_short', 'n_transcripts',
          'mean_f', 'f_eq_1', 'f_ge_0.90', 'f_0.10_to_0.90', 'f_lt_0.10', 'f_eq_0']


# ── Load all TSVs ──────────────────────────────────────────────────────────────
data     = {}
all_tids = set()
for fname, script in TSV_TO_SCRIPT.items():
    d = {}
    with open(RESULTS_DIR / fname) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            tid = row['transcript_id']
            d[tid] = float(row['uniqueness_factor'])
            all_tids.add(tid)
    data[script] = d

n_total = len(all_tids)
print(f"Total transcripts (union across all TSVs): {n_total}\n")

# ── Compute stats and write TSV ────────────────────────────────────────────────
out_path = RESULTS_DIR / 'results_summary.tsv'
with open(out_path, 'w', newline='') as fh:
    writer = csv.writer(fh, delimiter='\t')
    writer.writerow(FIELDS)

    for script in SCRIPT_ORDER:
        d    = data[script]
        fs   = list(d.values())
        n    = len(fs)
        short = n_total - n

        f_eq_1        = sum(1 for f in fs if f == 1.0)
        f_ge_90       = sum(1 for f in fs if f >= 0.90)
        f_lt_10       = sum(1 for f in fs if f <  0.10)
        f_eq_0        = sum(1 for f in fs if f == 0.0)
        f_mid         = n - f_eq_1 - f_lt_10          # 0.10 ≤ f < 1.0
        mean_f        = round(sum(fs) / n, 4) if n else 0

        ref, rsrc, mode, L, intron = METADATA[script]
        writer.writerow([script, ref, rsrc, mode, L, intron,
                         n_total, short, n,
                         mean_f, f_eq_1, f_ge_90, f_mid, f_lt_10, f_eq_0])

        check = short + f_eq_1 + f_mid + f_lt_10
        flag  = '' if check == n_total else f'  *** CHECK FAILS ({check} != {n_total})'
        print(f"  {script}  too_short={short:3d}  n={n}  mean={mean_f:.4f}  "
              f"f=1:{f_eq_1:>6}  f≥.90:{f_ge_90:>6}  "
              f"f.10-.90:{f_mid:>5}  f<.10:{f_lt_10:>5}  f=0:{f_eq_0:>5}{flag}")

print(f"\nSaved → {out_path}")
