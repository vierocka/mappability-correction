#!/usr/bin/env python3
"""transcript -> gene table, and one representative transcript per gene.

The representative is the transcript carrying --tag (MANE_Select for human GENCODE/Ensembl,
Ensembl_canonical for most other Ensembl species). Genes with no tagged transcript fall back to
the longest CDS, then the longest transcript, so the workflow runs on annotations that carry no
canonical tag at all.

usage: tx2gene.py --gtf FILE --out FILE --representatives FILE [--tag MANE_Select]
"""
import argparse
import gzip
import re
from collections import defaultdict

OPEN = lambda p: gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def attr(s, key):
    m = re.search(key + r'\s+"([^"]*)"', s)
    return m.group(1) if m else ""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gtf", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--representatives", required=True)
    p.add_argument("--tag", default="MANE_Select")
    a = p.parse_args()

    tx = {}
    exon_len = defaultdict(int)
    cds_len = defaultdict(int)
    with OPEN(a.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            feat = c[2]
            if feat not in ("transcript", "exon", "CDS"):
                continue
            tid = attr(c[8], "transcript_id")
            if not tid:
                continue
            if feat == "transcript":
                tx[tid] = {
                    "transcript_id": tid,
                    "gene_id": attr(c[8], "gene_id"),
                    "gene_name": attr(c[8], "gene_name") or attr(c[8], "gene_id"),
                    "biotype": (attr(c[8], "transcript_type") or attr(c[8], "transcript_biotype")
                                or attr(c[8], "gene_type") or attr(c[8], "gene_biotype")),
                    "tags": ";".join(re.findall(r'tag\s+"([^"]*)"', c[8])),
                }
            elif feat == "exon":
                exon_len[tid] += int(c[4]) - int(c[3]) + 1
            else:
                cds_len[tid] += int(c[4]) - int(c[3]) + 1

    if not tx:
        raise SystemExit("no transcript records found in the GTF")

    with open(a.out, "w") as fh:
        fh.write("transcript_id\tgene_id\tbiotype\tgene_name\ttranscript_len\tcds_len\ttags\n")
        for tid, t in tx.items():
            fh.write(f"{tid}\t{t['gene_id']}\t{t['biotype']}\t{t['gene_name']}\t"
                     f"{exon_len[tid]}\t{cds_len[tid]}\t{t['tags']}\n")

    by_gene = defaultdict(list)
    for tid, t in tx.items():
        by_gene[t["gene_id"]].append(tid)

    reps, n_tagged = [], 0
    for gid, tids in by_gene.items():
        tagged = [t for t in tids if a.tag and a.tag in tx[t]["tags"]]
        if tagged:
            n_tagged += 1
            reps.append(sorted(tagged, key=lambda t: -exon_len[t])[0])
        else:
            reps.append(sorted(tids, key=lambda t: (-cds_len[t], -exon_len[t], t))[0])

    with open(a.representatives, "w") as fh:
        fh.write("\n".join(sorted(reps)) + "\n")

    print(f"{len(tx)} transcripts, {len(by_gene)} genes; "
          f"{n_tagged} represented by a '{a.tag}' transcript, "
          f"{len(by_gene) - n_tagged} by the longest CDS/transcript")


if __name__ == "__main__":
    main()
