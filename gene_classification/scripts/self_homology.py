#!/usr/bin/env python3
"""Genes whose own sequence repeats internally.

A transcript built from a tandem repeat (polyubiquitin UBC, the NBPF DUF1220 arrays, mucin
tandem domains) loses uniqueness to itself, not to another locus. No external competitor is
involved, so scaling the count by funique would be wrong: the reads are real and belong to this
gene, they simply cannot be placed within it.

Measured as the fraction of a transcript's k-mers that occur more than once in the same
transcript, plus the longest internally repeated stretch.

usage: self_homology.py --fasta F --out F [--k 31] [--min-frac 0.05] [--tx2gene F]
"""
import argparse
import gzip
from collections import Counter

OPEN = lambda p: gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def read_fasta(path):
    name, chunks = None, []
    with OPEN(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def longest_repeat_run(seq, k, dup):
    """Longest run of consecutive positions whose k-mer is internally duplicated."""
    best = run = 0
    for i in range(len(seq) - k + 1):
        if seq[i:i + k] in dup:
            run += 1
            best = max(best, run)
        else:
            run = 0
    return best + k - 1 if best else 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--tx2gene")
    p.add_argument("--k", type=int, default=31)
    p.add_argument("--min-frac", type=float, default=0.05)
    a = p.parse_args()

    t2g = {}
    if a.tx2gene:
        with open(a.tx2gene) as fh:
            next(fh, None)
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 2:
                    t2g[c[0]] = (c[1], c[3] if len(c) > 3 else c[1])

    rows = []
    for tid, seq in read_fasta(a.fasta):
        seq = seq.upper()
        if len(seq) < a.k * 2:
            continue
        kmers = Counter(seq[i:i + a.k] for i in range(len(seq) - a.k + 1))
        dup = {km for km, n in kmers.items() if n > 1}
        if not dup:
            continue
        n_dup = sum(kmers[km] for km in dup)
        frac = n_dup / sum(kmers.values())
        if frac < a.min_frac:
            continue
        gid, gname = t2g.get(tid, (tid, tid))
        rows.append([gid, gname, tid, len(seq), round(frac, 4),
                     longest_repeat_run(seq, a.k, dup)])

    rows.sort(key=lambda r: -r[5])
    with open(a.out, "w") as fh:
        fh.write("gene_id\tgene_name\ttranscript_id\ttranscript_len\tfrac_repeated_kmers\tlongest_repeat_bp\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    print(f"{len(rows)} transcripts with >= {100*a.min_frac:.0f} % internally repeated {a.k}-mers")
    for r in rows[:10]:
        print(f"  {r[1]:<14} {r[5]:>7} bp repeat, {100*r[4]:.0f} % of k-mers")


if __name__ == "__main__":
    main()
