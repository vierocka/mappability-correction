#!/usr/bin/env python3
"""
07_sparse_celltypes/analyze.py
Task: Cell types at risk from marker signal loss in sparse scRNA-seq.

Context: Cell-type annotation in scRNA-seq relies on marker gene sets.  In sparse
10x Genomics data, a typical cell detects 20-40% of a cell type's markers.
If mappability artifacts affect some of those markers (uf_3p250 < 0.8 → biased or
zero signal), the detected subset shrinks further.  For cell types with few markers
(< 50), the combined effect can leave fewer than 5 reliable markers — making
confident annotation fragile.

This script:
  1. Loads CellMarker2 and PanglaoDB impact summaries
  2. Filters cell types with ≤ 200 markers in database (small marker sets)
  3. For each cell type, computes:
       expected_detected    = n_markers * detection_rate  (20% default)
       expected_uf_bad      = n_markers * frac_affected
       expected_reliable    = expected_detected - expected_detected * frac_affected
       risk_level: CRITICAL if expected_reliable < 2
                   HIGH     if expected_reliable < 5
                   MEDIUM   if expected_reliable < 10
  4. Lists the specific impaired markers per cell type

Inputs (independent):
  ../sc3UTRs/results/celltype_impact_summary.tsv
  ../sc3UTRs/results/marker_flags.tsv
  ../sc3UTRs/results/marker_3prime_impact.tsv    (used for per-gene details)

Outputs:
  sparse_celltype_risk.tsv      per cell type: risk level + expected reliable markers
  sparse_celltype_markers.tsv   for high/critical risk cell types: marker detail
  risk_summary.tsv              counts by risk level and database
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

HERE   = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from lib_check import check_library_tsv
ROOT   = HERE.parent.parent.parent
MAPDIR = ROOT / 'mappability'

CELLTYPE_IMPACT = MAPDIR / 'sc3UTRs' / 'results' / 'celltype_impact_summary.tsv'
MARKER_FLAGS    = MAPDIR / 'sc3UTRs' / 'results' / 'marker_flags.tsv'
MARKER_IMPACT   = MAPDIR / 'sc3UTRs' / 'results' / 'marker_3prime_impact.tsv'

# scRNA-seq detection parameters
DETECTION_RATE   = 0.20    # 20% of a cell type's markers detected per cell on average
# UF thresholds for "impaired" and "zero signal"
UF_IMPAIRED  = 0.80        # uf_3p250 < 0.8 → biased or reduced signal
UF_ZERO      = 0.10        # uf_3p250 < 0.1 → near-zero, unreliable

# Risk thresholds (expected reliable markers in detected set)
RISK_CRITICAL = 2
RISK_HIGH     = 5
RISK_MEDIUM   = 10

# Cell type size threshold
MAX_MARKERS = 200


def load_celltype_impact(path):
    """Returns list of rows from celltype_impact_summary.tsv."""
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            rows.append(row)
    return rows


def load_marker_flags(path):
    """Returns {gene_symbol: row} with worst UF metrics across all runs."""
    flags = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            flags[row['gene_symbol']] = row
    return flags


def load_marker_impact_v5g(path, label='v5g_MANE_RNA_PE_100'):
    """
    Load per-gene UF from the canonical run (v5g).
    Returns {gene_symbol: {uf_3p250: float, uf_overall: float, phenomenon: str}}
    """
    result = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row.get('label') != label:
                continue
            gsym = row.get('gene_symbol', '')
            if gsym:
                result[gsym] = {
                    'uf_3p250':   float(row.get('uf_3p250', 0) or 0),
                    'uf_overall': float(row.get('uf_overall', 0) or 0),
                    'phenomenon': row.get('phenomenon', ''),
                }
    return result


def risk_label(expected_reliable):
    if expected_reliable < RISK_CRITICAL:
        return 'CRITICAL'
    elif expected_reliable < RISK_HIGH:
        return 'HIGH'
    elif expected_reliable < RISK_MEDIUM:
        return 'MEDIUM'
    else:
        return 'LOW'


def main():
    # ── Library consistency check ─────────────────────────────────────────────
    # celltype_impact_summary.tsv aggregates marker impact across ALL runs
    # (multiple read lengths: 75/100/150/200 bp, SE/PE).
    # marker_flags.tsv stores worst_uf_3p250 = worst across all those runs.
    # marker_3prime_impact.tsv has per-run rows — we filter to label=v5g (PE 100 bp).
    #
    # IMPORTANT: the risk analysis uses worst_uf values (cross-run minimum).
    # This is conservative (takes the worst case regardless of read length).
    # UF is a fraction [0,1] — aggregating worst across lengths is valid as
    # an upper-bound risk estimate, but note that 75 bp SE gives lower UF than
    # 150 bp PE for the same gene (read length effect, not biology).
    # The 20% detection rate and frac_impaired are both ratios → normalisation OK.
    check_library_tsv(MARKER_IMPACT, label_col='label')

    print('Loading celltype impact summary …')
    impact_rows = load_celltype_impact(CELLTYPE_IMPACT)
    print(f'  {len(impact_rows):,} cell-type × database entries')

    print('Loading marker flags (worst-case per gene) …')
    marker_flags = load_marker_flags(MARKER_FLAGS)

    print('Loading per-gene UF from canonical run (v5g MANE RNA PE 100) …')
    gene_uf = load_marker_impact_v5g(MARKER_IMPACT)
    print(f'  {len(gene_uf):,} marker genes with UF data')

    # ── Filter small-marker cell types and compute risk ───────────────────────
    risk_rows = []
    detail_markers = []

    for row in impact_rows:
        n_markers = int(row.get('n_markers', 0))
        if n_markers == 0 or n_markers > MAX_MARKERS:
            continue

        db        = row.get('database', '')
        cell_type = row.get('cell_type', '')
        n_in_uf   = int(row.get('n_markers_in_uf', 0))
        n_affected = int(row.get('n_markers_affected', 0))

        # n_markers_affected in celltype_impact_summary used the overall 'affected' flag
        # which was: delta_3p500 < -0.10 OR uf_3p500 < 0.50
        # We recompute using uf_3p250 < 0.8 (tighter threshold, more relevant for 10x)

        affected_genes = [g.strip() for g in row.get('affected_markers', '').split(';') if g.strip()]

        # Recount using uf_3p250 threshold
        n_impaired  = 0   # uf_3p250 < 0.8
        n_zero      = 0   # uf_3p250 < 0.1
        impaired_genes = []
        zero_genes     = []
        for gsym in affected_genes:
            uf_data = gene_uf.get(gsym, {})
            uf3 = uf_data.get('uf_3p250', 1.0)   # default 1.0 if not found
            if uf3 < UF_ZERO:
                n_zero += 1
                zero_genes.append(gsym)
                n_impaired += 1
                impaired_genes.append(gsym)
            elif uf3 < UF_IMPAIRED:
                n_impaired += 1
                impaired_genes.append(gsym)

        # Fall back to n_in_uf if affected_genes list is empty
        if n_in_uf == 0:
            frac_impaired = 0.0
        else:
            frac_impaired = n_impaired / n_in_uf

        # Expected scenario: cell detects DETECTION_RATE fraction of markers
        expected_detected = n_markers * DETECTION_RATE
        # Of those detected, some fraction will have bad signal
        expected_bad       = expected_detected * frac_impaired
        expected_reliable  = expected_detected - expected_bad
        # Worst case: impaired + zero signal
        expected_zero_detected = n_markers * DETECTION_RATE * (n_zero / max(n_in_uf, 1))

        risk = risk_label(expected_reliable)

        risk_rows.append({
            'database':              db,
            'cell_type':             cell_type,
            'n_markers_db':          n_markers,
            'n_markers_in_uf_data':  n_in_uf,
            'n_impaired_uf3p250lt0.8': n_impaired,
            'n_zero_uf3p250lt0.1':   n_zero,
            'frac_impaired':         round(frac_impaired, 4),
            'expected_detected_20pct': round(expected_detected, 1),
            'expected_reliable':     round(expected_reliable, 1),
            'expected_zero_detected': round(expected_zero_detected, 1),
            'risk_level':            risk,
            'impaired_markers':      ';'.join(impaired_genes[:30]),   # cap for readability
            'zero_signal_markers':   ';'.join(zero_genes[:30]),
        })

        if risk in ('CRITICAL', 'HIGH'):
            for gsym in impaired_genes:
                uf_data = gene_uf.get(gsym, {})
                detail_markers.append({
                    'database':    db,
                    'cell_type':   cell_type,
                    'risk_level':  risk,
                    'gene_symbol': gsym,
                    'uf_3p250':    round(uf_data.get('uf_3p250', 1.0), 4),
                    'uf_overall':  round(uf_data.get('uf_overall', 1.0), 4),
                    'phenomenon':  uf_data.get('phenomenon', ''),
                    'flag_detail': marker_flags.get(gsym, {}).get('worst_uf_3p500', ''),
                })

    # Sort by risk then by expected_reliable
    risk_order = {'CRITICAL': 0, 'HIGH': 1, 'MEDIUM': 2, 'LOW': 3}
    risk_rows.sort(key=lambda r: (risk_order.get(r['risk_level'], 9), r['expected_reliable']))

    risk_fields = ['database', 'cell_type', 'risk_level',
                   'n_markers_db', 'n_markers_in_uf_data',
                   'n_impaired_uf3p250lt0.8', 'n_zero_uf3p250lt0.1',
                   'frac_impaired',
                   'expected_detected_20pct', 'expected_reliable', 'expected_zero_detected',
                   'impaired_markers', 'zero_signal_markers']
    with open(HERE / 'sparse_celltype_risk.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=risk_fields, delimiter='\t', lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(risk_rows)
    print(f'\nSaved sparse_celltype_risk.tsv  ({len(risk_rows)} cell types)')

    detail_fields = ['database', 'cell_type', 'risk_level', 'gene_symbol',
                     'uf_3p250', 'uf_overall', 'phenomenon', 'flag_detail']
    with open(HERE / 'sparse_celltype_markers.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=detail_fields, delimiter='\t', lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(detail_markers)
    print(f'Saved sparse_celltype_markers.tsv  ({len(detail_markers)} marker entries)')

    # ── Risk summary ──────────────────────────────────────────────────────────
    from collections import Counter
    risk_counts = Counter((r['database'], r['risk_level']) for r in risk_rows)
    dbs = sorted({r['database'] for r in risk_rows})
    levels = ['CRITICAL', 'HIGH', 'MEDIUM', 'LOW']
    summary_rows = []
    for db in dbs:
        for lvl in levels:
            n = risk_counts.get((db, lvl), 0)
            if n > 0:
                summary_rows.append({'database': db, 'risk_level': lvl, 'n_celltypes': n})

    with open(HERE / 'risk_summary.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['database', 'risk_level', 'n_celltypes'],
                           delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(summary_rows)
    print(f'Saved risk_summary.tsv')

    # ── Console ───────────────────────────────────────────────────────────────
    print('\n═══ Summary ═══')
    print(f'Cell types with ≤{MAX_MARKERS} markers: {len(risk_rows)}')
    for lvl in levels:
        n = sum(1 for r in risk_rows if r['risk_level'] == lvl)
        print(f'  {lvl:<10}: {n}')
    print(f'\nTop CRITICAL/HIGH cell types:')
    print(f'{"DB":<12} {"Cell type":<40} {"n_markers":>9} {"exp_reliable":>12} Risk')
    print('─' * 80)
    for r in risk_rows:
        if r['risk_level'] not in ('CRITICAL', 'HIGH'):
            continue
        ct = r['cell_type'][:38]
        print(f'{r["database"]:<12} {ct:<40} {r["n_markers_db"]:>9} '
              f'{r["expected_reliable"]:>12.1f} {r["risk_level"]}')


if __name__ == '__main__':
    main()
