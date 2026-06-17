#!/usr/bin/env python3
"""
shared_fraction.py — Quantify what fraction of each functional HLA transcript
                     is homologously shared with at least one pseudogene,
                     as a function of read-length window size.

Method:
  For each paralog group alignment (.aln.fa):
    1. Identify functional (F) and pseudogene (P) sequences.
    2. Slide a window of W bp (projected onto alignment columns) across the
       functional sequence, skipping columns that are a gap in that sequence.
    3. At each window, compute pairwise identity to every pseudogene in the group.
    4. A window is "shared" if max identity to any pseudogene >= THRESHOLD.
    5. Report: fraction of functional sequence windows that are shared.

Output:
  shared_fraction.tsv   — per functional gene, per window size, shared fraction
  shared_positions.tsv  — per-position shared flag (for plotting)
"""

import csv, re
from pathlib import Path
from collections import defaultdict

HERE      = Path(__file__).resolve().parent
THRESHOLDS = [0.90, 0.95]        # identity cutoffs
WINDOWS    = [75, 100, 150, 250]  # read-length window sizes (bp)

# Genes known to be pseudogenes (matched against sequence header tag)
# run_hla_mafft.py labels headers as:  gene (F|P) | ...
PSEUDO_TAG = '(P)'
FUNC_TAG   = '(F)'

GROUPS = [
    'classI_A_clade',
    'classI_BC_clade',
    'classII_DRB',
    'classII_DQ',
    'classII_DP',
]


# ── FASTA parsing ──────────────────────────────────────────────────────────────

def read_aligned_fasta(path):
    """Return list of (label, gene, is_pseudo, seq_with_gaps)."""
    records = []
    label = seq = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if label is not None:
                    records.append(label)
                label = line[1:]
                seq   = []
            else:
                seq.extend(line)
        if label is not None:
            records.append(label)   # flush — but seq not stored yet

    # Re-parse properly
    records = []
    label = seq_parts = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if label is not None:
                    gene      = label.split()[0]
                    is_pseudo = PSEUDO_TAG in label
                    records.append((label, gene, is_pseudo, ''.join(seq_parts)))
                label     = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if label is not None:
            gene      = label.split()[0]
            is_pseudo = PSEUDO_TAG in label
            records.append((label, gene, is_pseudo, ''.join(seq_parts)))
    return records


# ── Identity calculation ───────────────────────────────────────────────────────

def pairwise_identity(seq_a, seq_b):
    """
    Identity over aligned columns where BOTH sequences have a non-gap character.
    Returns (n_identical, n_comparable).
    """
    identical = comparable = 0
    for a, b in zip(seq_a, seq_b):
        if a == '-' or b == '-':
            continue
        comparable += 1
        if a.upper() == b.upper():
            identical += 1
    return identical, comparable


# ── Sliding window ─────────────────────────────────────────────────────────────

def alignment_columns_for_window(func_seq, start_base, window_bp):
    """
    Given the alignment sequence of the functional gene, find the range of
    alignment columns corresponding to [start_base, start_base + window_bp)
    in the ungapped functional sequence.
    Returns (col_start, col_end) or None if window exceeds sequence.
    """
    base_count = 0
    col_start  = None
    for col, ch in enumerate(func_seq):
        if ch == '-':
            continue
        if base_count == start_base:
            col_start = col
        base_count += 1
        if base_count == start_base + window_bp:
            return col_start, col + 1
    return None   # window extends past end of sequence


def shared_windows(func_seq, pseudo_seqs, window_bp, threshold):
    """
    Slide window_bp-base windows across func_seq (ungapped).
    Returns (n_shared, n_total) windows where at least one pseudogene
    has identity >= threshold over that window.
    """
    func_len  = sum(1 for c in func_seq if c != '-')
    n_total   = max(0, func_len - window_bp + 1)
    n_shared  = 0
    shared_at = []   # start positions (0-indexed in ungapped coords) that are shared

    for start in range(n_total):
        cols = alignment_columns_for_window(func_seq, start, window_bp)
        if cols is None:
            break
        col_s, col_e = cols
        f_window = func_seq[col_s:col_e]

        max_id = 0.0
        for p_seq in pseudo_seqs:
            p_window   = p_seq[col_s:col_e]
            ident, cmp = pairwise_identity(f_window, p_window)
            if cmp >= window_bp // 2:   # require at least half the window to be aligned
                max_id = max(max_id, ident / cmp)

        if max_id >= threshold:
            n_shared += 1
            shared_at.append(start)

    return n_shared, n_total, shared_at


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    rows_frac = []   # for shared_fraction.tsv
    rows_pos  = []   # for shared_positions.tsv (one row per position per condition)

    for group in GROUPS:
        aln_path = HERE / f'{group}.aln.fa'
        if not aln_path.exists():
            print(f'SKIP {group} — {aln_path.name} not found')
            continue

        records    = read_aligned_fasta(aln_path)
        func_recs  = [(lab, gene, seq) for lab, gene, is_p, seq in records if not is_p]
        pseudo_recs= [(lab, gene, seq) for lab, gene, is_p, seq in records if is_p]

        if not func_recs or not pseudo_recs:
            print(f'SKIP {group} — need both functional and pseudogene sequences')
            continue

        pseudo_seqs = [seq for _, _, seq in pseudo_recs]
        pseudo_genes = [gene for _, gene, _ in pseudo_recs]
        print(f'\n{group}')
        print(f'  Functional : {[g for _,g,_ in func_recs]}')
        print(f'  Pseudogenes: {pseudo_genes}')

        for _, func_gene, func_seq in func_recs:
            func_len = sum(1 for c in func_seq if c != '-')

            for W in WINDOWS:
                for T in THRESHOLDS:
                    n_shared, n_total, shared_at = shared_windows(
                        func_seq, pseudo_seqs, W, T)

                    frac = n_shared / n_total if n_total else 0.0
                    print(f'  {func_gene:<12}  W={W:>3} bp  threshold={T:.0%}  '
                          f'shared={n_shared}/{n_total}  ({100*frac:.1f}%)')

                    rows_frac.append({
                        'group':           group,
                        'functional_gene': func_gene,
                        'pseudogenes':     ';'.join(pseudo_genes),
                        'window_bp':       W,
                        'threshold_pct':   int(T * 100),
                        'n_windows_total': n_total,
                        'n_windows_shared':n_shared,
                        'frac_shared':     f'{frac:.4f}',
                        'func_len_bp':     func_len,
                    })

                    # Position track (only for W=100, T=90% to keep output manageable)
                    if W == 100 and T == 0.90:
                        shared_set = set(shared_at)
                        for pos in range(n_total):
                            rows_pos.append({
                                'group':           group,
                                'functional_gene': func_gene,
                                'position':        pos,
                                'shared_90pct':    int(pos in shared_set),
                            })

    # Write outputs
    frac_path = HERE / 'shared_fraction.tsv'
    with open(frac_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=[
            'group', 'functional_gene', 'pseudogenes',
            'window_bp', 'threshold_pct',
            'n_windows_total', 'n_windows_shared', 'frac_shared', 'func_len_bp'
        ], delimiter='\t')
        w.writeheader()
        w.writerows(rows_frac)
    print(f'\nSaved → {frac_path.name}')

    pos_path = HERE / 'shared_positions.tsv'
    with open(pos_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=[
            'group', 'functional_gene', 'position', 'shared_90pct'
        ], delimiter='\t')
        w.writeheader()
        w.writerows(rows_pos)
    print(f'Saved → {pos_path.name}')


if __name__ == '__main__':
    main()
