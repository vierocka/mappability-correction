#!/usr/bin/env python3
"""Are the most abundant transcripts the pseudogene-bearing ones?

Reads the pseudogene/uniqueness evidence table and ranks its transcripts by observed count in one
library. If inflation by pseudogene-derived reads matters for a net-sum repertoire, it has to show
up here: the net sum of a 9-mer is dominated by the most abundant transcripts that encode it.

Columns of resources/funique/low_uniqueness_MANE_transcripts.tsv used:
  UF_L75 / UF_L100              uniqueness factor at the library read length
  groups                        why uniqueness is low (pseudogenes / paralog / duplicate in MANE)
  pseudogene_copies_in_genome   pseudogene loci that compete for this transcript's reads
  named_pseudogenes             pseudogenes named after this gene in the annotation
  fcM_<group>_<run>             featureCounts -M --fraction count

Note the table holds only transcripts with low uniqueness (n = 1459), so it cannot be used to
compare against well-behaved genes; it answers "within the flagged set, where is the expression?".
Non-coding transcripts are removed by default, matching the pipeline: they carry no peptides, and in
total-RNA libraries the multicopy snRNAs (RNU1/RNU2) otherwise dominate the top of the ranking.

usage: inflation_evidence.py [-n 15] [--run fcM_JY_SRR364065] [--out DIR]
"""
import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import binomtest

HERE = Path(__file__).resolve().parent        # Collesano2026/pseudogene_problem/
ROOT = HERE.parent                            # Collesano2026/
TABLE = ROOT / "resources" / "funique" / "low_uniqueness_MANE_transcripts.tsv"

p = argparse.ArgumentParser()
p.add_argument("-n", type=int, default=15)
p.add_argument("--run", default="fcM_JY_SRR364065")
p.add_argument("--uf", default="UF_L100")
p.add_argument("--keep-non-coding", action="store_true",
               help="keep the 70 non-coding MANE transcripts; by default they are removed, as the "
                    "pipeline removes them before net sums (they carry no peptides)")
p.add_argument("--out", default=str(HERE))
a = p.parse_args()

NONCODING = ROOT / "resources" / "non_coding_MANE_transcripts_no_protein.tsv"

d = pd.read_csv(TABLE, sep="\t")
if not a.keep_non_coding:
    nc = set(pd.read_csv(NONCODING, sep="\t")["transcript_id"])
    n0 = len(d)
    d = d[~d.transcript_id.isin(nc)]
    print(f"removed {n0 - len(d)} non-coding transcripts (no protein -> no peptides)")
d["count"] = pd.to_numeric(d[a.run], errors="coerce")
d = d.dropna(subset=["count"])

d["has_pseudogene"] = d["groups"].fillna("").str.contains("pseudogene")
top = d.nlargest(a.n, "count")

cols = ["gene", a.uf, "pseudogene_copies_in_genome", "named_pseudogenes", "count", "groups"]
out = Path(a.out)
out.mkdir(parents=True, exist_ok=True)
dest = out / f"pseudogene_inflation_top_{a.run}.tsv"
top[cols].to_csv(dest, sep="\t", index=False)

print(top[cols].to_string(index=False, max_colwidth=44))
print(f"\nflagged transcripts in table: {len(d)}")
print(f"with pseudogene evidence:     {d.has_pseudogene.sum()} ({100*d.has_pseudogene.mean():.0f} %)")
print(f"of the top {a.n} by count:       {top.has_pseudogene.sum()} with pseudogene evidence")
print(f"median count, whole table:    {d['count'].median():,.0f}")
print(f"median count, top {a.n}:         {top['count'].median():,.0f}")
bt = binomtest(int(top.has_pseudogene.sum()), a.n, float(d.has_pseudogene.mean()), alternative="greater")
print(f"enrichment vs base rate:      binomial p = {bt.pvalue:.2g}")
print(f"\nwritten: {dest}")
