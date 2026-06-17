#!/usr/bin/env python3
"""
combo_summary.py — Summarise tool combinations detected across papers.

Reads:   tools_summary.tsv  (from mine_tools.py)

Output:
  combo_summary.tsv   — ranked list of tool combinations with counts
  combo_detail.tsv    — per-paper combo string for downstream joins
"""

import csv
from collections import Counter, defaultdict
from pathlib import Path

HERE   = Path(__file__).resolve().parent
IN     = HERE / 'tools_summary.tsv'
OUT_C  = HERE / 'combo_summary.tsv'
OUT_D  = HERE / 'combo_detail.tsv'

# Ordered within each category for reproducible combo strings
ALIGNERS = [
    'aln_STAR', 'aln_HISAT2', 'aln_TopHat', 'aln_bowtie', 'aln_subread',
    'aln_kallisto', 'aln_Salmon',
    'aln_CellRanger', 'aln_STARsolo', 'aln_Alevin', 'aln_kbpython',
    'aln_zUMIs', 'aln_Optimus',
]
COUNTERS = [
    'cnt_featureCounts', 'cnt_HTSeq', 'cnt_RSEM', 'cnt_bustools',
    'cnt_simpleaf', 'cnt_StringTie',
]
DOWNSTREAM = [
    'ds_DESeq2', 'ds_edgeR', 'ds_limma',
    'ds_Seurat', 'ds_Scanpy', 'ds_Monocle', 'ds_Harmony', 'ds_scVI',
]

# Human-readable short labels for combo string construction
SHORT = {
    'aln_STAR':          'STAR',
    'aln_HISAT2':        'HISAT2',
    'aln_TopHat':        'TopHat',
    'aln_bowtie':        'Bowtie',
    'aln_subread':       'Subread',
    'aln_kallisto':      'kallisto',
    'aln_Salmon':        'Salmon',
    'aln_CellRanger':    'CellRanger',
    'aln_STARsolo':      'STARsolo',
    'aln_Alevin':        'Alevin',
    'aln_kbpython':      'kb-python',
    'aln_zUMIs':         'zUMIs',
    'aln_Optimus':       'Optimus',
    'cnt_featureCounts': 'featureCounts',
    'cnt_HTSeq':         'HTSeq',
    'cnt_RSEM':          'RSEM',
    'cnt_bustools':      'bustools',
    'cnt_simpleaf':      'simpleaf',
    'cnt_StringTie':     'StringTie',
    'ds_DESeq2':         'DESeq2',
    'ds_edgeR':          'edgeR',
    'ds_limma':          'limma',
    'ds_Seurat':         'Seurat',
    'ds_Scanpy':         'Scanpy',
    'ds_Monocle':        'Monocle',
    'ds_Harmony':        'Harmony',
    'ds_scVI':           'scVI',
}

# Multimapper handling risk per aligner/counter
# True  = discards multimappers by default (signal loss in HLA / paralog-rich genes)
# False = keeps / redistributes multimappers (potential inflation)
# None  = depends on settings / not applicable
MULTI_RISK = {
    'aln_STAR':          'discard_default',   # --outSAMmultNmax 10, NH:i>1 reads in BAM
                                              # but featureCounts then drops NH>1 by default
    'aln_HISAT2':        'discard_default',   # same pattern with featureCounts
    'aln_TopHat':        'discard_default',
    'aln_bowtie':        'discard_default',
    'aln_subread':       'discard_default',
    'aln_kallisto':      'em_redistribute',   # EM splits ambiguous reads across targets
    'aln_Salmon':        'em_redistribute',
    'aln_CellRanger':    'discard_default',   # STARsolo default: NH=1 only counted
    'aln_STARsolo':      'discard_default',
    'aln_Alevin':        'em_redistribute',
    'aln_kbpython':      'em_redistribute',
    'aln_zUMIs':         'discard_default',
    'aln_Optimus':       'discard_default',
    'cnt_featureCounts': 'discard_default',   # -M flag required to count multimappers
    'cnt_HTSeq':         'discard_default',   # union mode: multimappers discarded
    'cnt_RSEM':          'em_redistribute',
    'cnt_bustools':      'em_redistribute',
    'cnt_simpleaf':      'em_redistribute',
    'cnt_StringTie':     'em_redistribute',
}


def true_keys(row, keys):
    return [k for k in keys if row.get(k, '').strip().lower() in ('true', '1', 'yes')]


def combo_string(aligners, counters, downstream):
    parts = (
        [SHORT[k] for k in aligners] +
        [SHORT[k] for k in counters] +
        [SHORT[k] for k in downstream]
    )
    return ' + '.join(parts) if parts else '(none detected)'


def classify_multi_risk(aligners, counters):
    all_tools = aligners + counters
    if not all_tools:
        return 'unknown'
    risk_types = {MULTI_RISK.get(t) for t in all_tools if t in MULTI_RISK}
    if 'em_redistribute' in risk_types and 'discard_default' in risk_types:
        return 'mixed'
    if 'em_redistribute' in risk_types:
        return 'em_redistribute'
    if 'discard_default' in risk_types:
        return 'discard_default'
    return 'unknown'


def main():
    rows = []
    with open(IN) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            rows.append(row)
    print(f'Loaded {len(rows)} papers from {IN.name}')

    combo_counts  = Counter()
    combo_meta    = {}   # combo_str → {n_scrna, n_bulk, years, multi_risk}
    detail_rows   = []
    year_combos   = defaultdict(Counter)

    for row in rows:
        aligners   = true_keys(row, ALIGNERS)
        counters   = true_keys(row, COUNTERS)
        downstream = true_keys(row, DOWNSTREAM)

        combo = combo_string(aligners, counters, downstream)
        multi = classify_multi_risk(aligners, counters)
        year  = row.get('year', '')

        combo_counts[combo] += 1
        year_combos[year][combo] += 1

        if combo not in combo_meta:
            combo_meta[combo] = {
                'n_scrna': 0, 'n_bulk': 0,
                'multi_risk': multi,
                'years': set(),
            }
        if row.get('likely_scrna', '').lower() in ('true', '1'):
            combo_meta[combo]['n_scrna'] += 1
        if row.get('likely_bulk', '').lower() in ('true', '1'):
            combo_meta[combo]['n_bulk'] += 1
        if year:
            combo_meta[combo]['years'].add(year)

        detail_rows.append({
            'doi':        row.get('doi', ''),
            'pmcid':      row.get('pmcid', ''),
            'year':       year,
            'journal':    row.get('journal', ''),
            'combo':      combo,
            'multi_risk': multi,
            'likely_scrna': row.get('likely_scrna', ''),
            'likely_bulk':  row.get('likely_bulk', ''),
            'bio_HLA':      row.get('bio_HLA', ''),
            'bio_MHC':      row.get('bio_MHC', ''),
            'bio_cancer':   row.get('bio_cancer', ''),
            'bio_immune_evasion': row.get('bio_immune_evasion', ''),
            'bio_DE_signature':   row.get('bio_DE_signature', ''),
            'bio_classify_patient': row.get('bio_classify_patient', ''),
            'bio_immunotherapy':    row.get('bio_immunotherapy', ''),
            'bio_pseudogene':       row.get('bio_pseudogene', ''),
            'github_links': row.get('github_links', ''),
        })

    n_total = len(rows)
    years   = sorted(year_combos)

    # ── Write combo_summary.tsv ───────────────────────────────────────────────
    fields = ['combo', 'n_papers', 'pct', 'n_scrna', 'n_bulk', 'multi_risk', 'year_range'] + years
    with open(OUT_C, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for combo, n in combo_counts.most_common():
            meta = combo_meta[combo]
            yr_range = (f"{min(meta['years'])}–{max(meta['years'])}"
                        if meta['years'] else '')
            out = {
                'combo':      combo,
                'n_papers':   n,
                'pct':        f'{100*n/n_total:.1f}',
                'n_scrna':    meta['n_scrna'],
                'n_bulk':     meta['n_bulk'],
                'multi_risk': meta['multi_risk'],
                'year_range': yr_range,
            }
            for yr in years:
                out[yr] = year_combos[yr].get(combo, 0)
            w.writerow(out)
    print(f'Saved {len(combo_counts)} unique combos → {OUT_C}')

    # ── Write combo_detail.tsv ────────────────────────────────────────────────
    det_fields = [
        'doi', 'pmcid', 'year', 'journal', 'combo', 'multi_risk',
        'likely_scrna', 'likely_bulk',
        'bio_HLA', 'bio_MHC', 'bio_cancer', 'bio_immune_evasion',
        'bio_DE_signature', 'bio_classify_patient', 'bio_immunotherapy',
        'bio_pseudogene', 'github_links',
    ]
    with open(OUT_D, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=det_fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(detail_rows)
    print(f'Saved {len(detail_rows)} rows → {OUT_D}')

    # ── Console summary ───────────────────────────────────────────────────────
    print(f'\nTop 20 tool combinations (n={n_total} papers total):')
    print(f'\n{"Combo":<55} {"n":>5} {"pct":>6}  {"risk":<18} {"scRNA":>5} {"bulk":>5}')
    print('─' * 100)
    for combo, n in combo_counts.most_common(20):
        pct  = f'{100*n/n_total:.1f}%'
        meta = combo_meta[combo]
        print(f'{combo:<55} {n:>5} {pct:>6}  {meta["multi_risk"]:<18} '
              f'{meta["n_scrna"]:>5} {meta["n_bulk"]:>5}')


if __name__ == '__main__':
    main()
