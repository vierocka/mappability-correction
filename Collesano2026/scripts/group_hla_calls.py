#!/usr/bin/env python3
"""Runs of one group (cell line / individual) -> one HLA genotype.
Takes the run with most HLA reads (sum of gene_reads); reports per-run concordance."""
import argparse

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--group", required=True)
ap.add_argument("--calls", nargs="+", required=True, help="<run>_hla_calls.tsv")
ap.add_argument("--out", required=True)
ap.add_argument("--concordance", required=True)
a = ap.parse_args()

t = pd.concat([pd.read_csv(f, sep="\t") for f in a.calls])
t["genotype"] = t.allele1 + "/" + t.allele2
best = t.groupby("sample").gene_reads.sum().idxmax()
final = t[t["sample"] == best].copy()
final.insert(0, "group", a.group)
final["source_run"] = best
final.to_csv(a.out, sep="\t", index=False)

w = t.pivot(index="gene", columns="sample", values="genotype")
w["final"] = final.set_index("gene").genotype
w["all_runs_agree"] = w.drop(columns="final").nunique(axis=1) == 1
w.insert(0, "group", a.group)
w.to_csv(a.concordance, sep="\t")
if not w.all_runs_agree.all():
    print(f"WARNING {a.group}: runs disagree for {list(w.index[~w.all_runs_agree])}; using {best}")
