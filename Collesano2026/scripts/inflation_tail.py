#!/usr/bin/env python3
"""
Does x funique act mainly on the extreme (inflated) right tail of 9-mer net sums?
  before = exclusions applied, no funique (TPM renormalised)
  after  = TPM_corrected (pipeline output)
  touched = peptide with >= 1 source transcript corrected (f < 1)
  n_src   = number of MANE transcripts encoding the peptide (multicopy)
Outputs: _tail_bins.tsv, _tail_summary.tsv, _hist.png, _shift.png
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--ninemers", required=True)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    per, f = [], None
    for p in a.samples:
        d = pd.read_csv(p, sep="\t").set_index("gene_id")
        rpk = (d.counts_corrected / d.f_unique) / (d.length / 1e3)
        per.append(pd.DataFrame({"before": rpk / rpk.sum() * 1e6, "after": d.TPM_corrected}))
        f = d.f_unique if f is None else np.fmin(f, d.f_unique)
    t = sum(per) / len(per)
    t["corrected"] = f < 1

    k = pd.read_csv(a.ninemers, sep="\t", usecols=["gene_id", "peptide_9mer", "count_in_gene"]).join(t, on="gene_id")
    k = k[k.before.notna()]
    g = k.groupby("peptide_9mer")
    pep = pd.DataFrame({
        "before": (k.before * k.count_in_gene).groupby(k.peptide_9mer).sum(),
        "after": (k.after * k.count_in_gene).groupby(k.peptide_9mer).sum(),
        "touched": g.corrected.any(),
        "n_src": g.size(),
    })
    pep = pep[pep.before > 0]
    pep["log2_ratio"] = np.log2(pep.after.clip(lower=1e-12) / pep.before)
    pct = pep.before.rank(pct=True)
    edges = [0, .5, .9, .99, .999, 1.0]
    labels = ["0-50%", "50-90%", "90-99%", "99-99.9%", "top 0.1%"]
    pep["bin"] = pd.cut(pct, edges, labels=labels, include_lowest=True)
    b = pep.groupby("bin", observed=True).agg(
        n=("before", "size"), pct_touched=("touched", "mean"), pct_multicopy=("n_src", lambda x: (x > 1).mean()),
        median_log2_ratio_touched=("log2_ratio", lambda x: x[pep.loc[x.index, "touched"]].median()),
        min_before=("before", "min"))
    b["pct_touched"] *= 100; b["pct_multicopy"] *= 100
    b.insert(0, "cell_line", a.cell)
    b.to_csv(f"{a.out_prefix}_tail_bins.tsv", sep="\t", float_format="%.4g")

    tt = pep[pep.touched]
    s = pd.Series({
        "cell_line": a.cell, "peptides": len(pep), "touched": len(tt),
        "pct_touched": 100 * len(tt) / len(pep),
        "pct_of_touched_in_top1pct": 100 * (pct[tt.index] > .99).mean(),
        "pct_of_touched_in_top10pct": 100 * (pct[tt.index] > .9).mean(),
        "median_log2_ratio_touched": tt.log2_ratio.median(),
        "median_log2_ratio_untouched": pep[~pep.touched].log2_ratio.median(),
        "netsum_mass_removed_pct": 100 * (1 - pep.after.sum() / pep.before.sum()),
        "max_before": pep.before.max(), "max_after": pep.after.max(),
        "top10_before": ",".join(pep.before.nlargest(10).index),
        "top10_after": ",".join(pep.after.nlargest(10).index),
    })
    s.to_frame("value").to_csv(f"{a.out_prefix}_tail_summary.tsv", sep="\t")

    bins = np.linspace(-3, np.log10(pep.before.max()) + .2, 80)
    # below 2:1 overall: a wider figure is scaled to the page width when rendered and its
    # axis labels become unreadable
    fig, ax = plt.subplots(1, 2, figsize=(9, 4.4))
    for col, c in [("before", "#999999"), ("after", "#1f77b4")]:
        ax[0].hist(np.log10(pep[col].clip(lower=1e-3)), bins=bins, histtype="step", color=c, label=col)
        ax[1].hist(np.log10(tt[col].clip(lower=1e-3)), bins=bins, histtype="step", color=c, label=col)
    for x, ttl in zip(ax, ["all 9-mers", f"touched 9-mers (n={len(tt):,})"]):
        x.set_yscale("log"); x.set_xlabel("log10 net-sum TPM"); x.set_title(f"{a.cell}: {ttl}", fontsize=9)
        x.legend(frameon=False, fontsize=8)
    ax[0].set_ylabel("peptides")
    plt.tight_layout(); plt.savefig(f"{a.out_prefix}_hist.png", dpi=150); plt.close()

    fig, ax = plt.subplots(figsize=(5, 3.6))
    hb = ax.hexbin(np.log10(tt.before), tt.log2_ratio, gridsize=60, bins="log", cmap="Blues", mincnt=1)
    ax.set_xlabel("log10 net-sum before"); ax.set_ylabel("log2(after / before)")
    ax.set_title(f"{a.cell}: touched 9-mers", fontsize=9); plt.colorbar(hb, label="peptides")
    plt.tight_layout(); plt.savefig(f"{a.out_prefix}_shift.png", dpi=150); plt.close()
    print(b.to_string()); print(s.drop(["top10_before", "top10_after"]).to_string())


if __name__ == "__main__":
    main()
