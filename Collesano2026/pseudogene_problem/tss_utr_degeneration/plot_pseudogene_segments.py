#!/usr/bin/env python3
"""Figure for pseudogene_segment_conservation.tsv. usage: plot_pseudogene_segments.py OUTDIR"""
import csv, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

OUT = Path(sys.argv[1])
rows = {r["pseudogene"]: r for r in csv.DictReader(open(OUT / "pseudogene_segment_conservation.tsv"), delimiter="\t")}
GROUPS = [("full copy", ["HLA-H", "HLA-J", "HLA-K"]),
          ("ends decaying", ["HLA-L", "HLA-W", "HLA-DPB2", "MICE"]),
          ("promoter lost", ["HLA-DRB6", "HLA-DPA2", "HLA-DQB3"])]
COLS = ["promoter_1kb", "UTR5", "CDS", "UTR3"]
order = [g for _, gs in GROUPS for g in gs]
M = np.array([[float(rows[g][c]) for c in COLS] for g in order])

ink, muted = "#1f1f1d", "#6b6b66"
cmap = LinearSegmentedColormap.from_list("blue", ["#f4f8fd", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
cmap.set_bad("#e6e5e1")
fig, ax = plt.subplots(figsize=(7.2, 6.0))
ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=0, vmax=1, aspect="auto")
for i in range(M.shape[0]):
    for j in range(4):
        v = M[i, j]
        ax.text(j, i, "too short" if np.isnan(v) else f"{v:.2f}", ha="center", va="center",
                fontsize=8 if np.isnan(v) else 10, color=muted if np.isnan(v) else ("white" if v >= 0.55 else ink))
# white gaps between groups
y = -0.5
for name, gs in GROUPS:
    ax.text(3.62, y + len(gs) / 2, name, rotation=270, va="center", ha="left", fontsize=10, color=ink, fontweight="bold")
    y += len(gs)
    if y < len(order) - 0.5: ax.axhline(y, color="white", lw=4)
ax.set_xticks(range(4), ["promoter\n(1 kb)", "5′UTR", "CDS", "3′UTR"], fontsize=10, color=ink)
ax.set_yticks(range(len(order)), [f"{g}  ← {rows[g]['parent']}" for g in order], fontsize=10, color=ink)
ax.tick_params(length=0); ax.xaxis.tick_top()
for s in ax.spines.values(): s.set_visible(False)
fig.suptitle("MHC pseudogenes: same protein code, different regulatory ends", x=0.02, ha="left",
             fontsize=13, fontweight="bold", color=ink)
fig.text(0.02, 0.915, "Share of each part of the functional gene still found in the pseudogene locus (0 = gone, 1 = identical).\n"
         "GENCODE v49, GRCh38, blastn. Grey: UTR shorter than 50 bp.", fontsize=8.5, color=muted, ha="left", va="top")
fig.subplots_adjust(left=0.27, right=0.90, top=0.79, bottom=0.03)
fig.savefig(OUT / "pseudogene_segment_conservation.png", dpi=220)
