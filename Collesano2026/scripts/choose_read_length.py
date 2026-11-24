#!/usr/bin/env python3
"""
Per-mate mean read length (parse_fastqc.py output) -> funique table length.
  73-77 bp  -> 75   (precomputed table)
  98-102 bp -> 100  (precomputed table)
  else      -> round(mean)  (table simulated by the pipeline, once per length)
Writes: '<L>\t<mean_len>\t<source>' (source = precomputed / simulate).
"""
import argparse

import pandas as pd

PRECOMPUTED = {75: (73, 77), 100: (98, 102)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", required=True, help="parse_fastqc.py TSV (rows = mates) of the reads that are aligned")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    m = float(pd.read_csv(a.qc, sep="\t").mean_length.mean())
    L, src = round(m), "simulate"
    for k, (lo, hi) in PRECOMPUTED.items():
        if lo <= m <= hi:
            L, src = k, "precomputed"
    open(a.out, "w").write(f"{L}\t{m:.2f}\t{src}\n")
    print(f"mean read length {m:.2f} bp -> funique L{L} ({src})")


if __name__ == "__main__":
    main()
