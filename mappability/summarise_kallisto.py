#!/usr/bin/env python3
"""
Summarise Kallisto pseudoalignment results across all read source × index combinations.

Output tables (written to results/):
  kallisto_run_summary.tsv   — per-run pseudoalignment rates from run_info.json
  kallisto_inflation.tsv     — per-transcript TPM: MANE index vs dedup CDS index,
                               for each read source × mode × read length

Inflation direction:
  ratio_mane_over_cds > 1  →  MANE index inflates (UTR k-mers rescue ambiguous reads)
  ratio_mane_over_cds < 1  →  CDS index inflates (CDS unique but MANE transcript multimaps)
  ratio = NaN              →  CDS has no match in id_conversion_table (XP_/YP_ or non-coding)

CDS reads → CDS self-map show only ~5–6% unique reads despite 100% pseudoalignment:
the coding sequences are highly redundant but MANE's UTRs mask this by providing
unique k-mers that absorb the multimapping ambiguity.
"""

import csv
import json
import re
from collections import defaultdict
from pathlib import Path

HERE        = Path(__file__).resolve().parent
KALLISTO    = HERE / "kallisto_results"
RESULTS     = HERE / "results"
ID_TABLE    = RESULTS / "id_conversion_table.tsv"

# ── Parse directory name into metadata ────────────────────────────────────────

def parse_dir(name):
    """Return (read_source, mode, read_len) from mappability_* directory name."""
    # dedup CDS reads
    if "dedup_cds_reads" in name:
        src = "CDS_reads"
    elif "exonflank" in name or "exon_flank" in name:
        src = "exon-flank"
    elif "MANE_RNA" in name or "MANE_exonflank" not in name and "MANE" in name:
        src = "MANE_RNA"
    else:
        src = "unknown"

    # override for exon-flank
    if "exonflank" in name:
        src = "exon-flank"

    m = re.search(r'_L(\d+)bp', name)
    read_len = m.group(1) if m else "100"

    mode = "PE" if "_PE" in name else "SE"

    return src, mode, read_len

# ── Build NP_ → ENST map from id_conversion_table ────────────────────────────

np_to_enst   = {}
enst_to_gene = {}
with open(ID_TABLE) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        enst = row['ensembl_transcript_id']
        np   = row['refseq_protein_id']
        sym  = row['gene_symbol']
        if np:
            np_to_enst[np] = enst
        enst_to_gene[enst] = sym

print(f"ID table: {len(np_to_enst):,} NP_ → ENST mappings, {len(enst_to_gene):,} ENST entries")

# ── Helper: parse NP_ from dedup CDS target_id ────────────────────────────────
# Format: lcl|NC_000001.11_cds_NP_001005484.2_1
#      or lcl|NC_000001.11_cds_XP_047292308.1_2

_np_re = re.compile(r'_cds_([A-Z]+_\d+\.\d+)_\d+$')

def extract_protein_id(target_id):
    m = _np_re.search(target_id)
    return m.group(1) if m else None

# ── Load abundance.tsv ────────────────────────────────────────────────────────

def load_abundance_mane(path):
    """Returns {enst: tpm}"""
    d = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            d[row['target_id']] = float(row['tpm'])
    return d

def load_abundance_cds_by_enst(path):
    """Returns {enst: tpm_sum} — aggregates all CDS entries per NP_, maps to ENST."""
    counts = defaultdict(float)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            np_id = extract_protein_id(row['target_id'])
            if np_id and np_id in np_to_enst:
                enst = np_to_enst[np_id]
                counts[enst] += float(row['tpm'])
    return dict(counts)

# ── Collect all run directories ───────────────────────────────────────────────

run_dirs = sorted(p for p in KALLISTO.iterdir() if p.is_dir())

# ── Output 1: run-level summary ───────────────────────────────────────────────

run_summary_path = RESULTS / "kallisto_run_summary.tsv"
run_fields = ['read_dir', 'read_source', 'mode', 'read_len_bp', 'index',
              'n_processed', 'n_pseudoaligned', 'p_pseudoaligned',
              'n_unique', 'p_unique']

run_rows = []
for rdir in run_dirs:
    src, mode, rlen = parse_dir(rdir.name)
    for idx_name in ('kallisto_vs_MANE', 'kallisto_vs_dedup_cds'):
        ri_path = rdir / idx_name / "run_info.json"
        if not ri_path.exists():
            continue
        ri = json.loads(ri_path.read_text())
        run_rows.append({
            'read_dir':        rdir.name,
            'read_source':     src,
            'mode':            mode,
            'read_len_bp':     rlen,
            'index':           'MANE' if 'MANE' in idx_name else 'dedup_CDS',
            'n_processed':     ri['n_processed'],
            'n_pseudoaligned': ri['n_pseudoaligned'],
            'p_pseudoaligned': ri['p_pseudoaligned'],
            'n_unique':        ri['n_unique'],
            'p_unique':        ri['p_unique'],
        })

with open(run_summary_path, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=run_fields, delimiter='\t')
    w.writeheader()
    w.writerows(run_rows)

print(f"\nRun summary ({len(run_rows)} rows) → {run_summary_path}")
print(f"\n{'Read source':<20} {'Mode':>4} {'Len':>5}  "
      f"{'MANE p_align':>13} {'MANE p_uniq':>11}  "
      f"{'CDS p_align':>12} {'CDS p_uniq':>10}")
print("─" * 80)
for rdir in run_dirs:
    src, mode, rlen = parse_dir(rdir.name)
    mane_ri = rdir / 'kallisto_vs_MANE'     / 'run_info.json'
    cds_ri  = rdir / 'kallisto_vs_dedup_cds'/ 'run_info.json'
    if not mane_ri.exists() or not cds_ri.exists():
        continue
    m = json.loads(mane_ri.read_text())
    c = json.loads(cds_ri.read_text())
    print(f"{src:<20} {mode:>4} {rlen:>5}  "
          f"{m['p_pseudoaligned']:>12.1f}% {m['p_unique']:>10.1f}%  "
          f"{c['p_pseudoaligned']:>11.1f}% {c['p_unique']:>9.1f}%")

# ── Output 2: per-transcript inflation ───────────────────────────────────────

infl_path  = RESULTS / "kallisto_inflation.tsv"
infl_fields = ['transcript_id', 'gene_symbol',
               'read_source', 'mode', 'read_len_bp',
               'tpm_mane', 'tpm_cds', 'ratio_mane_over_cds', 'direction']

infl_rows = []
for rdir in run_dirs:
    src, mode, rlen = parse_dir(rdir.name)

    mane_ab = rdir / 'kallisto_vs_MANE'      / 'abundance.tsv'
    cds_ab  = rdir / 'kallisto_vs_dedup_cds' / 'abundance.tsv'
    if not mane_ab.exists() or not cds_ab.exists():
        continue

    tpm_mane = load_abundance_mane(mane_ab)
    tpm_cds  = load_abundance_cds_by_enst(cds_ab)

    all_enst = set(tpm_mane) | set(tpm_cds)
    for enst in sorted(all_enst):
        tm = tpm_mane.get(enst, 0.0)
        tc = tpm_cds.get(enst, 0.0)

        if tm == 0.0 and tc == 0.0:
            continue   # no signal in either — skip

        if tc > 0:
            ratio = round(tm / tc, 4)
            direction = 'MANE_higher' if tm > tc else ('CDS_higher' if tc > tm else 'equal')
        else:
            ratio     = None
            direction = 'CDS_zero'

        infl_rows.append({
            'transcript_id':      enst,
            'gene_symbol':        enst_to_gene.get(enst, ''),
            'read_source':        src,
            'mode':               mode,
            'read_len_bp':        rlen,
            'tpm_mane':           round(tm, 4),
            'tpm_cds':            round(tc, 4),
            'ratio_mane_over_cds': ratio if ratio is not None else 'NA',
            'direction':          direction,
        })

with open(infl_path, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=infl_fields, delimiter='\t')
    w.writeheader()
    w.writerows(infl_rows)

print(f"\nInflation table ({len(infl_rows):,} rows) → {infl_path}")

# ── Per-run direction summary ─────────────────────────────────────────────────

from collections import Counter

print(f"\n{'Read source':<20} {'Mode':>4} {'Len':>5}  "
      f"{'n_transcripts':>14}  {'MANE_higher':>12}  {'CDS_higher':>11}  {'CDS_zero':>9}")
print("─" * 85)

# Group by (src, mode, rlen)
groups = defaultdict(list)
for row in infl_rows:
    key = (row['read_source'], row['mode'], row['read_len_bp'])
    groups[key].append(row['direction'])

for key in sorted(groups):
    src, mode, rlen = key
    cnt = Counter(groups[key])
    n   = sum(cnt.values())
    print(f"{src:<20} {mode:>4} {rlen:>5}  "
          f"{n:>14,}  "
          f"{cnt['MANE_higher']:>12,}  "
          f"{cnt['CDS_higher']:>11,}  "
          f"{cnt['CDS_zero']:>9,}")

print(f"\nSaved:\n  {run_summary_path}\n  {infl_path}")
