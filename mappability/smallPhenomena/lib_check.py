"""
lib_check.py — Read-length and library-size consistency checks.

Import and call check_library() at the top of every analyze.py before any
statistics.  The check reports what read lengths and conditions are present
in the input files and warns if mixing is detected.

Normalization notes
-------------------
  UF (uniqueness_factor / uf_overall / uf_3p250 etc.)
      Fraction [0, 1].  Already normalised by position count within each run.
      Safe to compare across different read lengths — but be aware that longer
      reads produce higher UF in general (more sequence context → more unique
      alignments).  Comparisons across read lengths are directional, not metric.

  TPM (kallisto tpm column)
      Transcripts Per Million.  Normalised for transcript length and library
      size within each run.  Safe to compare across samples of the SAME read
      length and reference.  Cross-length TPM comparison is risky (different
      effective lengths).

  est_counts (kallisto est_counts column)
      Raw EM-estimated counts.  NOT normalised for library size.
      Do NOT compare est_counts across runs without dividing by total counts
      (or converting to TPM first).

Usage
-----
    from lib_check import check_library

    # For a single known run:
    check_library(
        run_label   = 'v5g_MANE_RNA_PE_100',
        read_len_bp = 100,
        mode        = 'PE',
        n_reads     = None,    # fill if known from run_info.json
    )

    # For a TSV with a 'label' or 'read_len' column (checks all values present):
    check_library_tsv(tsv_path, label_col='label')
"""

import csv
import re
import sys
from collections import Counter
from pathlib import Path


def check_library(run_label, read_len_bp, mode='PE', n_reads=None,
                  metrics_are_pct=True):
    """
    Print a brief library summary and flag any concerns.

    Parameters
    ----------
    run_label   : str   Human-readable run name (e.g. 'v5g_MANE_RNA_PE_100')
    read_len_bp : int   Read length in bp
    mode        : str   'PE' or 'SE'
    n_reads     : int or None   Total reads in the library (from run_info.json)
    metrics_are_pct : bool   If True, prints confirmation that metrics are %-based
    """
    print(f'\n[lib_check] Run: {run_label}')
    print(f'            Read length: {read_len_bp} bp  |  Mode: {mode}')
    if n_reads is not None:
        print(f'            Library size: {n_reads:,} reads')
    if metrics_are_pct:
        print(f'            Metrics: UF and uf_3pW are fractions [0,1]. '
              f'TPM is per-million. Safe for cross-sample % comparison.')
    print()


def check_library_tsv(tsv_path, label_col='label', len_col=None):
    """
    Scan a results TSV for read-length / condition metadata.
    Reports all unique values found in label_col and warns if multiple
    read lengths are mixed.

    Parameters
    ----------
    tsv_path  : Path or str
    label_col : str   Column containing condition/run labels (e.g. 'label')
    len_col   : str or None   Column with numeric read length if present
    """
    path = Path(tsv_path)
    if not path.exists():
        print(f'[lib_check] WARNING: {path.name} not found — skipping check.')
        return

    labels   = Counter()
    lengths  = Counter()
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            if label_col in row:
                labels[row[label_col]] += 1
            if len_col and len_col in row:
                try:
                    lengths[int(row[len_col])] += 1
                except ValueError:
                    pass

    # Try to infer read length from label names (e.g. "v5g_MANE_RNA_PE_100")
    inferred = set()
    for lbl in labels:
        m = re.search(r'_(\d{2,3})(?:bp)?(?:_|$)', lbl)
        if m:
            inferred.add(int(m.group(1)))

    all_lens = lengths.keys() | inferred

    print(f'[lib_check] File: {path.name}')
    print(f'            Unique conditions: {len(labels)}')
    if labels:
        for lbl, n in sorted(labels.items())[:5]:
            print(f'              {lbl}  ({n} rows)')
        if len(labels) > 5:
            print(f'              … and {len(labels)-5} more')

    if len(all_lens) > 1:
        print(f'  WARNING  Multiple read lengths detected: {sorted(all_lens)}')
        print(f'           UF comparisons across lengths are directional, not metric.')
        print(f'           TPM comparisons across lengths require length normalisation.')
        print(f'           est_counts MUST NOT be compared across lengths without')
        print(f'           converting to TPM or dividing by total library counts.')
    elif all_lens:
        print(f'            Read length: {next(iter(all_lens))} bp — single condition, OK.')
    print()


def warn_if_mixing_runs(file_paths_and_labels):
    """
    Given a list of (file_path, read_len_bp) tuples, warn if multiple read
    lengths are about to be merged into a single statistics table.

    Example
    -------
        warn_if_mixing_runs([
            (STAR_UF, 100),
            (KAL_ABUND, 100),
        ])
    """
    lens = Counter(rlen for _, rlen in file_paths_and_labels)
    if len(lens) > 1:
        print('[lib_check] WARNING: mixing files from different read lengths:')
        for path, rlen in file_paths_and_labels:
            print(f'  {rlen} bp  {Path(path).name}')
        print('  UF comparisons: directional only.  TPM: normalise first.  '
              'est_counts: do not compare.')
    else:
        rlen = next(iter(lens))
        print(f'[lib_check] All input files: {rlen} bp — consistent, OK.')
    print()
