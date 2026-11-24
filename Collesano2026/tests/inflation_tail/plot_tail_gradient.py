#!/usr/bin/env python3
"""The correction's effect against expression rank, faceted by cell line.

Two quantities per net-sum bin, both read from *_tail_bins.tsv:
  pct_touched    — share of 9-mers in the bin whose net sum changes under x funique
  pct_multicopy  — share of 9-mers in the bin that come from more than one source gene
                   (a multicopy gene family)

Both rise monotonically with expression, which is the claim the figure has to carry: the
correction is confined to the most abundant peptides, so a global rank statistic over 9.4 M
9-mers cannot see it.

usage: plot_tail_gradient.py [--outdir DIR]
"""
import argparse
import glob
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ORDER = ["0-50%", "50-90%", "90-99%", "99-99.9%", "top 0.1%"]
# one hue per series, assigned in fixed order and never cycled
C_TOUCHED, C_MULTI = "#3b6fb6", "#c2622d"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default=".")
    a = p.parse_args()
    here = pathlib.Path(__file__).resolve().parent
    out = pathlib.Path(a.outdir)

    srcs = sorted(glob.glob(str(here / "*_tail_bins.tsv")))
    if not srcs:
        raise SystemExit("no *_tail_bins.tsv next to this script")
    d = pd.concat([pd.read_csv(f, sep="\t") for f in srcs])
    d["bin"] = pd.Categorical(d["bin"], categories=ORDER, ordered=True)
    d = d.sort_values(["cell_line", "bin"])
    cells = sorted(d.cell_line.unique())

    # keep the panels side by side for comparison, but below 2:1 overall: a wider figure is
    # scaled to the page width when rendered and its labels become unreadable
    fig, axes = plt.subplots(1, len(cells), figsize=(3.2 * len(cells), 5.0),
                             sharey=True, constrained_layout=True)
    axes = [axes] if len(cells) == 1 else list(axes)
    x = range(len(ORDER))

    for ax, cell in zip(axes, cells):
        g = d[d.cell_line == cell].set_index("bin").reindex(ORDER)
        ax.plot(x, g.pct_touched, "-o", lw=2, ms=6, color=C_TOUCHED, label="changed by × funique")
        ax.plot(x, g.pct_multicopy, "-s", lw=2, ms=6, color=C_MULTI, label="from a multicopy family")
        # label the end point of each series: the number the figure exists to show
        for s, c in ((g.pct_touched, C_TOUCHED), (g.pct_multicopy, C_MULTI)):
            ax.annotate(f"{s.iloc[-1]:.0f} %", (x[-1], s.iloc[-1]), textcoords="offset points",
                        xytext=(-4, 9), ha="right", color=c, fontsize=9, fontweight="bold")
        ax.set_title(cell, fontsize=11)
        ax.set_xticks(list(x))
        ax.set_xticklabels(ORDER, rotation=45, ha="right", fontsize=8)
        ax.set_xlabel("net-sum percentile bin", fontsize=9)
        ax.grid(axis="y", lw=.4, alpha=.35)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    axes[0].set_ylabel("% of 9-mers in the bin", fontsize=9)
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle("The mappability correction is confined to the most abundant peptides",
                 fontsize=12, y=1.06)

    for ext in ("png", "pdf"):
        fig.savefig(out / f"tail_gradient.{ext}", dpi=200, bbox_inches="tight")
    print(f"wrote {out / 'tail_gradient.png'} and .pdf from {len(srcs)} cell lines")

    piv = d.pivot_table(index="bin", columns="cell_line", values=["pct_touched", "pct_multicopy"],
                        observed=False).reindex(ORDER)
    print()
    print(piv.round(2).to_string())


if __name__ == "__main__":
    main()
