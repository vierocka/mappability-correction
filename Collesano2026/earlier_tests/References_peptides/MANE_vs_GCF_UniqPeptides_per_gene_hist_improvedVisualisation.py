#!/usr/bin/env python3
"""
plot_peptide_distributions.py
==============================
Two-panel figure:
  Top    — step (outline) histograms, log/log, MANE vs GCF
  Bottom — per-bin difference (GCF − MANE gene count), log x / symlog y

Input  : mane_9mers.tsv  /  gcf_9mers.tsv   (output of peptides_per_gene.py)
Output : peptide_distribution.png  (300 dpi)

Usage  : python plot_peptide_distributions.py
"""

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── Config ────────────────────────────────────────────────────────────────────

MANE_TSV = "mane_9mers.tsv"
GCF_TSV  = "gcf_9mers.tsv"
OUT_PNG  = "peptide_distribution.png"
BINS     = 50

C_MANE = "#0072B2"   # blue  (Wong 2011)
C_GCF  = "#E69F00"   # amber (Wong 2011)
C_POS  = "#009E73"   # green  — GCF > MANE
C_NEG  = "#D55E00"   # red-orange — GCF < MANE

# ── Load data ─────────────────────────────────────────────────────────────────

def load_counts(path):
    counts = []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                counts.append(int(parts[2]))
    return np.array(counts, dtype=np.int64)

print("Loading MANE …", file=sys.stderr)
mane = load_counts(MANE_TSV)
print("Loading GCF …",  file=sys.stderr)
gcf  = load_counts(GCF_TSV)

# ── Shared log-spaced bin edges ───────────────────────────────────────────────

combined_min = max(1, min(mane.min(), gcf.min()))
combined_max = max(mane.max(), gcf.max())
bin_edges = np.logspace(np.log10(combined_min),
                        np.log10(combined_max),
                        BINS + 1)

# Pre-compute histogram counts for both (needed for the difference panel)
mane_counts, _ = np.histogram(mane, bins=bin_edges)
gcf_counts,  _ = np.histogram(gcf,  bins=bin_edges)
diff_counts    = gcf_counts.astype(np.int64) - mane_counts.astype(np.int64)
bin_centres    = np.sqrt(bin_edges[:-1] * bin_edges[1:])   # geometric midpoint

# ── Figure layout: 2 rows ─────────────────────────────────────────────────────

fig, (ax_top, ax_bot) = plt.subplots(
    2, 1,
    figsize=(10, 8),
    gridspec_kw={'height_ratios': [3, 1.4], 'hspace': 0.08},
    sharex=True
)

# ════════════════════════════════════════════════════════════════════════════
# TOP PANEL — step histograms
# ════════════════════════════════════════════════════════════════════════════

ax_top.hist(mane, bins=bin_edges, histtype='step',
            color=C_MANE, linewidth=2.0,
            label=f"MANE canonical  (n={len(mane):,} genes)", zorder=4)

ax_top.hist(gcf, bins=bin_edges, histtype='step',
            color=C_GCF, linewidth=2.0, linestyle='--',
            label=f"GCF all isoforms (n={len(gcf):,} genes)", zorder=3)

# Filled area for GCF only where it exceeds MANE (highlight surplus)
gcf_surplus = np.where(gcf_counts > mane_counts, gcf_counts, np.nan)
ax_top.stairs(gcf_surplus, bin_edges,
              fill=True, color=C_GCF, alpha=0.18, zorder=2,
              label="GCF surplus over MANE")

ax_top.set_yscale('log')

# Median lines
for vals, color, label in [(mane, C_MANE, f"MANE median: {int(np.median(mane)):,}"),
                            (gcf,  C_GCF,  f"GCF median: {int(np.median(gcf)):,}")]:
    ax_top.axvline(np.median(vals), color=color, linewidth=1.4,
                   linestyle=':', zorder=5, label=label)

ax_top.set_ylabel("Number of genes (log₁₀)", fontsize=11)
ax_top.yaxis.set_major_formatter(
    ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
ax_top.yaxis.set_minor_formatter(ticker.NullFormatter())
ax_top.spines['top'].set_visible(False)
ax_top.spines['right'].set_visible(False)
ax_top.grid(axis='both', which='major', linewidth=0.5, alpha=0.35, zorder=0)
ax_top.grid(axis='both', which='minor', linewidth=0.25, alpha=0.15, zorder=0)
ax_top.legend(fontsize=9.5, framealpha=0.92, loc='upper left')

# Stats box
def stats_str(label, arr):
    return (f"{label}\n"
            f"  median : {int(np.median(arr)):>8,}\n"
            f"  mean   : {int(np.mean(arr)):>8,}\n"
            f"  min    : {int(arr.min()):>8,}\n"
            f"  max    : {int(arr.max()):>8,}\n"
            f"  genes  : {len(arr):>8,}")

ax_top.text(0.99, 0.98,
            stats_str("MANE", mane) + "\n\n" + stats_str("GCF", gcf),
            transform=ax_top.transAxes, fontsize=7.5,
            va='top', ha='right', fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                      edgecolor='#cccccc', alpha=0.92))

# ════════════════════════════════════════════════════════════════════════════
# BOTTOM PANEL — per-bin difference GCF − MANE
# ════════════════════════════════════════════════════════════════════════════

# Bar colours: green where GCF > MANE, red-orange where GCF < MANE
bar_colors = [C_POS if d >= 0 else C_NEG for d in diff_counts]

ax_bot.bar(bin_centres, diff_counts,
           width=np.diff(bin_edges),   # actual bin widths (log-scaled)
           color=bar_colors, alpha=0.80,
           align='center', zorder=3)

ax_bot.axhline(0, color='#444444', linewidth=0.9, zorder=4)

# symlog so zero is visible and both positive/negative are log-scaled
max_abs = max(abs(diff_counts.max()), abs(diff_counts.min()), 1)
ax_bot.set_yscale('symlog', linthresh=1)
ax_bot.set_ylim(-max_abs * 2, max_abs * 2)

ax_bot.set_ylabel("GCF − MANE\ngenes per bin", fontsize=10)
ax_bot.set_xlabel("Unique 9-mer peptides per gene (log₁₀)", fontsize=11)

ax_bot.yaxis.set_major_formatter(
    ticker.FuncFormatter(lambda x, _: f"{int(x):+,}"))
ax_bot.xaxis.set_major_formatter(
    ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
ax_bot.xaxis.set_minor_formatter(ticker.NullFormatter())

ax_bot.spines['top'].set_visible(False)
ax_bot.spines['right'].set_visible(False)
ax_bot.grid(axis='both', which='major', linewidth=0.5, alpha=0.35, zorder=0)

# Colour legend for difference panel
from matplotlib.patches import Patch
ax_bot.legend(handles=[Patch(color=C_POS, alpha=0.8, label="GCF > MANE"),
                        Patch(color=C_NEG, alpha=0.8, label="GCF < MANE")],
              fontsize=9, framealpha=0.9, loc='upper right')

# Shared x log scale (set on bottom; top shares via sharex)
ax_bot.set_xscale('log')

plt.savefig(OUT_PNG, dpi=300, bbox_inches='tight')
print(f"Saved → {OUT_PNG}", file=sys.stderr)
