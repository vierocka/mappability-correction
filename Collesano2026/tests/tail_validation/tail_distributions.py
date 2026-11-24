#!/usr/bin/env python3
"""Net-sum distributions before and after correction, and where MS ligands sit in them.

Three questions the tail test depends on:
  1. does the correction change the SHAPE of the net-sum distribution, or only its scale?
  2. where in the distribution do MS-detected ligands actually fall?
  3. in the tail, do the before and after rankings differ at all?

Writes a quantile table, an MS-position table and one figure. Reads only; nothing is written
outside this folder.

usage: tail_distributions.py --data-dir DIR [--outdir .]
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tail_validation import MEAN_COUNTS, MS_FILES, load_ms, net_sums

QS = [0.5, 0.9, 0.99, 0.999, 0.9999, 1.0]
C_BEFORE, C_AFTER, C_MS = "#3b6fb6", "#c2622d", "#4a4a4a"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data-dir", required=True)
    p.add_argument("--cells", nargs="+", default=["JY", "HCC1937", "HCT116"])
    p.add_argument("--outdir", default=".")
    a = p.parse_args()
    data = pathlib.Path(a.data_dir).expanduser()
    out = pathlib.Path(a.outdir)

    qrows, mrows = [], []
    fig, axes = plt.subplots(2, len(a.cells), figsize=(4.0 * len(a.cells), 6.4),
                             constrained_layout=True)
    if len(a.cells) == 1:
        axes = axes.reshape(2, 1)

    for j, cell in enumerate(a.cells):
        ms = load_ms(data, cell)
        pep, _ = net_sums(data, cell)
        is_ms = pep.index.isin(ms)

        for q in QS:
            qrows.append({"cell_line": cell, "quantile": q,
                          "before": float(pep.before.quantile(q)),
                          "after": float(pep.after.quantile(q)),
                          "ratio_after_before": float(pep.after.quantile(q) / pep.before.quantile(q))})

        pct = pep.after.rank(pct=True) * 100
        mrows.append({"cell_line": cell, "n_ms": int(is_ms.sum()),
                      "ms_median_percentile": float(pct[is_ms].median()),
                      "ms_q25_percentile": float(pct[is_ms].quantile(.25)),
                      "ms_q75_percentile": float(pct[is_ms].quantile(.75)),
                      "ms_pct_above_p99": float(100 * (pct[is_ms] > 99).mean()),
                      "ms_pct_above_p999": float(100 * (pct[is_ms] > 99.9).mean()),
                      "ms_median_netsum": float(pep.after[is_ms].median()),
                      "all_median_netsum": float(pep.after.median())})

        # top: full distribution, before vs after, with the MS ligands' positions
        ax = axes[0, j]
        lo = max(pep.before[pep.before > 0].min(), 1e-3)
        bins = np.logspace(np.log10(lo), np.log10(pep.before.max()) + .2, 70)
        ax.hist(pep.before.clip(lower=lo), bins=bins, histtype="step", lw=2,
                color=C_BEFORE, label="before")
        ax.hist(pep.after.clip(lower=lo), bins=bins, histtype="step", lw=2,
                color=C_AFTER, label="after x funique")
        ax.hist(pep.after[is_ms].clip(lower=lo), bins=bins, histtype="stepfilled", lw=0,
                color=C_MS, alpha=.55, label=f"MS ligands (n={is_ms.sum()})",
                weights=np.full(is_ms.sum(), len(pep) / max(is_ms.sum(), 1) / 200))
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(cell, fontsize=11)
        ax.set_xlabel("net sum (TPM)", fontsize=9)
        if j == 0:
            ax.set_ylabel("9-mers per bin (log)", fontsize=9)
            ax.legend(frameon=False, fontsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

        # bottom: the tail only — after vs before, per-peptide
        ax = axes[1, j]
        tail = pep.nlargest(int(len(pep) * .001), "before")
        ax.scatter(tail.before, tail.after, s=4, alpha=.25, color=C_BEFORE, lw=0,
                   label="top 0.1 % (before)")
        if is_ms.sum():
            t_ms = tail[tail.index.isin(ms)]
            ax.scatter(t_ms.before, t_ms.after, s=26, color=C_MS, lw=0,
                       label=f"MS ligands in tail (n={len(t_ms)})")
        lim = [tail.before.min(), tail.before.max()]
        ax.plot(lim, lim, "--", lw=1, color="#999999", label="unchanged")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("net sum before", fontsize=9)
        if j == 0:
            ax.set_ylabel("net sum after", fontsize=9)
            ax.legend(frameon=False, fontsize=8, loc="upper left")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    fig.suptitle("Net-sum distributions and the top 0.1 % tail, before and after x funique",
                 fontsize=12)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"tail_distributions.{ext}", dpi=190, bbox_inches="tight")

    q = pd.DataFrame(qrows); m = pd.DataFrame(mrows)
    q.to_csv(out / "netsum_quantiles.tsv", sep="\t", index=False)
    m.to_csv(out / "ms_position.tsv", sep="\t", index=False)
    print("quantiles of the net-sum distribution:")
    print(q.to_string(index=False))
    print("\nwhere MS ligands sit (percentile of the corrected net sum):")
    print(m.to_string(index=False))
    print(f"\nwrote {out/'tail_distributions.png'}, netsum_quantiles.tsv, ms_position.tsv")


if __name__ == "__main__":
    main()
