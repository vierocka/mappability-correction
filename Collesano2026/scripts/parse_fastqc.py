#!/usr/bin/env python3
"""
FastQC fastqc_data.txt (text) -> one QC row per file.
  mean_length   from 'Sequence Length Distribution' (ranges -> midpoint)
  mode_length   most frequent length
  mean_phred    from 'Per sequence quality scores' (read-weighted)
  phred_ok      mean_phred >= --min-phred (default 30)
STAR Log.final.out (optional) adds 'Average input read length' (both mates).
"""
import argparse
import re
import sys


def module(lines, name):
    out, on = [], False
    for ln in lines:
        if ln.startswith(f">>{name}"):
            on = True
            continue
        if on and ln.startswith(">>END_MODULE"):
            break
        if on and not ln.startswith("#"):
            out.append(ln.rstrip("\n").split("\t"))
    return out


def mid(x):
    a = [float(v) for v in x.split("-")]
    return sum(a) / len(a)


def parse(path):
    L = open(path).read().splitlines()
    basic = dict(r[:2] for r in module(L, "Basic Statistics") if len(r) >= 2)
    ld = [(mid(r[0]), float(r[1])) for r in module(L, "Sequence Length Distribution")]
    n = sum(c for _, c in ld)
    q = [(float(r[0]), float(r[1])) for r in module(L, "Per sequence quality scores")]
    nq = sum(c for _, c in q)
    return {
        "file": basic.get("Filename", path),
        "encoding": basic.get("Encoding", ""),
        "total_sequences": int(basic.get("Total Sequences", 0)),
        "pct_gc": basic.get("%GC", ""),
        "mean_length": round(sum(l * c for l, c in ld) / n, 2) if n else None,
        "mode_length": max(ld, key=lambda t: t[1])[0] if ld else None,
        "mean_phred": round(sum(v * c for v, c in q) / nq, 2) if nq else None,
    }


def star_len(path):
    for ln in open(path):
        if "Average input read length" in ln:
            return float(ln.split("|")[1])
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastqc", nargs="+", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--stage", default="raw", help="raw / trimmed")
    ap.add_argument("--min-phred", type=float, default=30)
    ap.add_argument("--star-log", default=None)
    ap.add_argument("--out", default="-")
    a = ap.parse_args()
    rows = [parse(p) for p in a.fastqc]
    sl = star_len(a.star_log) if a.star_log else None
    cols = ["sample", "stage", "file", "encoding", "total_sequences", "pct_gc", "mean_length",
            "mode_length", "mean_phred", "phred_ok", "star_avg_input_len_pair"]
    out = sys.stdout if a.out == "-" else open(a.out, "w")
    out.write("\t".join(cols) + "\n")
    for r in rows:
        r.update(sample=a.sample, stage=a.stage, star_avg_input_len_pair=sl,
                 phred_ok=(r["mean_phred"] or 0) >= a.min_phred)
        out.write("\t".join(str(r[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
