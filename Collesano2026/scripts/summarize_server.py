#!/usr/bin/env python3
"""Server pipeline summaries: HLA types, read QC (length, Phred, STAR), corrected TPM matrix, overview plot."""
import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--samples", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--min-phred", type=float, default=30)
a = ap.parse_args()
O = Path(a.outdir); S = O / "summary"; S.mkdir(exist_ok=True)
sm = pd.read_csv(a.samples, sep="\t", dtype=str).fillna("")

# HLA
hla = pd.concat([pd.read_csv(O / f"groups/{g}/hla/{g}_hla_calls.tsv", sep="\t") for g in sm.group.unique()])
hla[["group", "gene", "allele1", "allele2", "decision", "optitype", "arcasHLA", "support_final", "gene_reads",
     "cellosaurus_best", "cellosaurus_match", "source_run"]].to_csv(S / "hla_types.tsv", sep="\t", index=False)


def star(path):
    d = {}
    for ln in open(path):
        if "|" in ln:
            k, v = ln.split("|"); d[k.strip()] = v.strip()
    return d


rows = []
for r in sm.itertuples():
    for st in ("raw", "trimmed"):
        q = pd.read_csv(O / f"runs/{r.run}/qc/{r.run}_qc_{st}.tsv", sep="\t")
        rows.append({"run": r.run, "group": r.group, "stage": st, "reads": int(q.total_sequences.iloc[0]),
                     "mean_length": round(q.mean_length.mean(), 2), "mean_phred": round(q.mean_phred.mean(), 2),
                     "phred_ok": bool(q.mean_phred.mean() >= a.min_phred)})
    L, m, src = open(O / f"runs/{r.run}/qc/{r.run}_read_length.txt").read().split()
    s = star(O / f"runs/{r.run}/align/{r.run}_{r.group}_MANEpersUTR_Log.final.out")
    rows[-1].update(funique_table=f"L{L} ({src})")
    rows[-2].update(funique_table=f"L{L} ({src})",
                    star_avg_input_len=s.get("Average input read length"),
                    star_unique_pct=s.get("Uniquely mapped reads %"),
                    star_multi_pct=s.get("% of reads mapped to multiple loci"))
qc = pd.DataFrame(rows); qc.to_csv(S / "qc_reads.tsv", sep="\t", index=False)

# corrected TPM matrix (transcript x run) + group/condition means
cols = {}
for r in sm.itertuples():
    d = pd.read_csv(O / f"correction/{r.group}/corrected_counts/{r.run}_mane_fc_corrected.tsv", sep="\t", index_col=0)
    cols[f"{r.run}_TPM_corrected"] = d.TPM_corrected
for (g, c), _ in sm.groupby(["group", "condition"], sort=False):
    d = pd.read_csv(O / f"correction/{g}/corrected_counts/{g}_corrected_mean_{c}_mane_fc_corrected.tsv",
                    sep="\t", index_col=0)
    cols[f"{g}_{c}_mean_TPM_corrected"] = d.mean_TPM_corrected
pd.DataFrame(cols).to_csv(S / "corrected_TPM_matrix.tsv", sep="\t", float_format="%.6f")

# overview plot
hv = pd.read_csv(O / "correction/HLA_corrected_mean_values.tsv", sep="\t")
raw = qc[qc.stage == "raw"]
fig, ax = plt.subplots(1, 3, figsize=(14, 3.8))
ax[0].bar(raw.run, raw.mean_phred, color=["#1f77b4" if ok else "#d62728" for ok in raw.phred_ok])
ax[0].axhline(a.min_phred, ls="--", c="k", lw=.8); ax[0].set_ylabel("mean Phred (raw)")
ax[1].bar(raw.run, raw.mean_length, color="#7f7f7f"); ax[1].set_ylabel("mean read length per mate (bp)")
for x in (ax[0], ax[1]):
    x.tick_params(axis="x", rotation=90, labelsize=7)
GENE = {"ENST00000376809": "HLA-A", "ENST00000412585": "HLA-B", "ENST00000376228": "HLA-C",
        "ENST00000354258": "TAP1", "ENST00000374897": "TAP2"}
key = hv.iloc[:, 0].str.split(":")
hv["group"] = key.str[0].str.split("/").str[0]
hv["gene"] = key.str[1].str.split(".").str[0].map(GENE).fillna("other")
p = hv.pivot_table(index="group", columns="gene", values="mean_TPM_corrected", aggfunc="sum")
p.plot.bar(ax=ax[2], width=.8); ax[2].set_ylabel("corrected TPM (allele sum)"); ax[2].legend(fontsize=7, frameon=False)
ax[2].tick_params(axis="x", rotation=0)
ax[0].set_title("read quality", fontsize=9); ax[1].set_title("read length → funique table", fontsize=9)
ax[2].set_title("HLA / TAP corrected TPM", fontsize=9)
plt.tight_layout(); plt.savefig(S / "overview.png", dpi=150)
print(hla.to_string(index=False))
