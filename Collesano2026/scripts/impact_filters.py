#!/usr/bin/env python3
"""
Impact of each processing step on the 9-mer net-sum table, scored against MS ligands.
n_source_repertoire_9mers: 9-mers enumerated from the expressed personal MANE proteome
(net sum > 0); not a binding prediction.
Variants (transcript level, then 9-mer net sum = sum TPM x count_in_gene):
  raw            TPM from fcM counts, all transcripts
  excl           exclusions removed (+alt copies folded), TPM renormalised, no funique
  excl_uf        + x funique  (= TPM_corrected, pipeline output)
  *_min1         same, transcripts with TPM < 1 dropped
Per variant: n peptides, MS recall, ROC-AUC (MS vs rest), Spearman vs excl_uf, top-1% overlap.
"""
import argparse
import re

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def auc(score, label):
    r = pd.Series(score).rank().values
    n1 = label.sum(); n0 = len(label) - n1
    return (r[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--samples", nargs="+", required=True, help="<SRR>_mane_fc_corrected.tsv of one condition")
    ap.add_argument("--ninemers", required=True, help="9-mer table (gene_id, gene_symbol, peptide_9mer, count_in_gene)")
    ap.add_argument("--table-s1", required=True)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    per = []
    for p in a.samples:                                   # variants per sample, then replicate mean
        d = pd.read_csv(p, sep="\t").set_index("gene_id")
        cnt = (d.counts_corrected / d.f_unique)           # undo x funique, keep alt-copy fold
        rpk = cnt / (d.length / 1e3)
        per.append(pd.DataFrame({"raw": d.TPM, "excl": rpk / rpk.sum() * 1e6, "excl_uf": d.TPM_corrected}))
    tpm = sum(per) / len(per)
    for v in ["raw", "excl", "excl_uf"]:
        tpm[f"{v}_min1"] = tpm[v].where(tpm[v] >= 1)

    k = pd.read_csv(a.ninemers, sep="\t", usecols=["gene_id", "peptide_9mer", "count_in_gene"])
    k = k.join(tpm, on="gene_id")
    pep = {}
    for v in tpm.columns:
        w = k[v] * k.count_in_gene
        pep[v] = w.groupby(k.peptide_9mer).sum(min_count=1)
    pep = pd.DataFrame(pep)

    s1 = pd.read_excel(a.table_s1, sheet_name="Peptides")
    cols = [c for c in s1.columns if re.fullmatch(rf"Intensity {re.escape(a.cell)}-\d+", c)]
    ms = set(s1.loc[s1[cols].fillna(0).gt(0).any(axis=1) & (s1.Length == 9), "Sequence"])

    ref = pep.excl_uf.fillna(0)
    top_ref = set(ref.nlargest(int(len(ref) * .01)).index)
    rows = []
    for v in pep.columns:
        x = pep[v].fillna(0)
        pos = x[x > 0]
        lab = pos.index.isin(ms)
        rows.append({"cell_line": a.cell, "variant": v, "n_source_repertoire_9mers": len(pos),
                     "ms_9mers": len(ms), "recall_pct": round(100 * lab.sum() / len(ms), 2),
                     "auc": round(auc(pos.values, lab), 4),
                     "spearman_vs_excl_uf": round(spearmanr(x, ref).correlation, 4),
                     "top1pct_overlap_pct": round(100 * len(set(x.nlargest(len(top_ref)).index) & top_ref) / len(top_ref), 2),
                     "sum_netsum": x.sum()})
    out = pd.DataFrame(rows)
    out.to_csv(f"{a.out_prefix}_filter_impact.tsv", sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
