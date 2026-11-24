#!/usr/bin/env python3
"""MANE (or personal MANE) transcript FASTA -> featureCounts SAF: one feature per transcript, 1..length.
(cds_mane.saf: the name is historical; features span the whole transcript, UTRs included.)"""
import sys

src, dst = sys.argv[1], sys.argv[2]
rows, tid, n = [], None, 0
for ln in open(src):
    if ln.startswith(">"):
        if tid:
            rows.append((tid, n))
        tid, n = ln[1:].split()[0], 0
    else:
        n += len(ln.strip())
if tid:
    rows.append((tid, n))
with open(dst, "w") as f:
    f.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    for t, L in rows:
        f.write(f"{t}\t{t}\t1\t{L}\t+\n")
