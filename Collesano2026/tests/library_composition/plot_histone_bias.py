#!/usr/bin/env python3
"""The chemistry marker that survives non-coding removal and rescaling.

Replication-dependent histones are the one protein-coding family without a poly(A) tail, so their
share of coding TPM reports library chemistry rather than biology. Plots that share per library
from `library_composition.tsv`; the per-gene table is used only for the spread.

usage: plot_histone_bias.py [--indir .] [--outdir .]
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

C = {"poly(A)": "#3b6fb6", "Ribo-Zero": "#c2622d"}   # as in tail_distributions.py


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--indir", default=".")
    p.add_argument("--outdir", default=".")
    a = p.parse_args()
    d = pathlib.Path(a.indir)
    out = pathlib.Path(a.outdir)

    lib = pd.read_csv(d / "library_composition.tsv", sep="\t")
    gene = pd.read_csv(d / "histone_tpm_by_library.tsv", sep="\t").set_index("gene")

    lib = lib.sort_values(["chemistry", "RD_histone_pct_of_PC_TPM"])
    y = lib.RD_histone_pct_of_PC_TPM.values
    lo, hi = y[lib.chemistry == "poly(A)"], y[lib.chemistry == "Ribo-Zero"]

    fig, ax = plt.subplots(figsize=(6.4, 3.8), constrained_layout=True)
    x = np.arange(len(lib))
    ax.bar(x, y, width=.62, color=[C[c] for c in lib.chemistry], zorder=3)
    for xi, v in zip(x, y):
        ax.text(xi, v * 1.12, f"{v:.3f}", ha="center", va="bottom", fontsize=9)

    ax.set_yscale("log")
    ax.set_ylim(0.02, 5)
    ax.set_xticks(x)
    ax.set_xticklabels(lib.library, fontsize=9)
    ax.set_ylabel("RD histones, % of protein-coding TPM", fontsize=9)
    ax.set_title(f"{hi.min()/lo.max():.0f}–{hi.max()/lo.min():.0f}x more histone signal in the "
                 "Ribo-Zero library,\nafter non-coding removal and rescaling", fontsize=10)
    ax.grid(axis="y", lw=.5, color="#dddddd", zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=C[c]) for c in ("poly(A)", "Ribo-Zero")]
    ax.legend(handles, ["poly(A)", "Ribo-Zero"], frameon=False, fontsize=9, loc="upper left")

    for ext in ("png", "pdf"):
        fig.savefig(out / f"histone_bias.{ext}", dpi=190, bbox_inches="tight")

    n_det = int((gene > 0).any(axis=1).sum())
    print(f"{len(gene)} RD histone genes ({n_det} detected in at least one library)")
    print(lib[["library", "chemistry", "RD_histone_pct_of_PC_TPM"]].to_string(index=False))
    print(f"fold difference: {hi.min()/lo.max():.1f}-{hi.max()/lo.min():.1f}x")
    print(f"wrote {out/'histone_bias.png'}")


if __name__ == "__main__":
    main()
