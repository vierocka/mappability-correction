#!/usr/bin/env python3
"""
plot_shared_fraction.py — 3-panel figure: shared sequence fraction for HLA-A, DRB1, B
                           with 5'UTR / CDS / 3'UTR annotation track per panel.

Reads:   shared_positions.tsv  (from shared_fraction.py, W=100 bp, threshold=90%)
Output:  hla_shared_fraction.pdf  +  .png
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from pathlib import Path

HERE = Path(__file__).resolve().parent

GENES = ['HLA-A', 'HLA-DRB1', 'HLA-B']

PSEUDOGENES = {
    'HLA-A':    'pseudogenes: H, J, L, P, V',
    'HLA-DRB1': 'pseudogenes: DRB2, DRB6, DRB7, DRB8, DRB9',
    'HLA-B':    'pseudogenes: K, N, T, U, W',
}

# UTR/CDS boundaries in transcript coordinates (0-indexed)
# Derived from GENCODE v49 GTF for MANE Select transcripts:
#   HLA-A   ENST00000376809  + strand  1535 bp
#   HLA-DRB1 ENST00000360004  - strand  1223 bp
#   HLA-B   ENST00000412585  - strand  1536 bp
BOUNDARIES = {
    'HLA-A':    {'utr5_end': 21,  'cds_start': 22,  'utr3_start': 1117, 'tx_len': 1535},
    'HLA-DRB1': {'utr5_end': 105, 'cds_start': 106, 'utr3_start': 904,  'tx_len': 1223},
    'HLA-B':    {'utr5_end': 20,  'cds_start': 21,  'utr3_start': 1107, 'tx_len': 1536},
}

SMOOTH = 80   # rolling mean window (bp)

# Colours
COL_SHARED = '#b03a2e'
COL_LINE   = '#7b241c'
COL_UNIQUE = '#d5d8dc'
COL_50     = '#566573'
COL_UTR    = '#7fb3d3'   # 5' and 3' UTR — muted blue
COL_CDS    = '#2e4053'   # CDS — dark slate
COL_BOUND  = '#1a252f'   # boundary tick lines


def smooth(arr, win):
    return (pd.Series(arr)
              .rolling(win, center=True, min_periods=5)
              .mean()
              .values)


def draw_annotation(ax, b, xmax):
    """
    Draw a gene-structure annotation bar inside ax using data coordinates.
    The bar occupies y = 1.04 – 1.13; vertical boundaries drawn from 0 to bar top.
    xmax: last window-start position (= tx_len - 100).
    """
    bar_y0, bar_y1 = 1.04, 1.13
    bh = bar_y1 - bar_y0

    cds_s  = b['cds_start']
    utr3_s = b['utr3_start']

    # 5'UTR block
    ax.add_patch(mpatches.FancyArrow(
        0, (bar_y0 + bar_y1) / 2, cds_s, 0,
        width=bh * 0.55, head_width=bh * 0.55, head_length=0,
        color=COL_UTR, lw=0, zorder=6,
    ))
    ax.add_patch(mpatches.Rectangle(
        (0, bar_y0), cds_s, bh, color=COL_UTR, lw=0, zorder=6,
    ))

    # CDS block
    ax.add_patch(mpatches.Rectangle(
        (cds_s, bar_y0), utr3_s - cds_s, bh,
        color=COL_CDS, lw=0, zorder=6,
    ))

    # 3'UTR block
    ax.add_patch(mpatches.Rectangle(
        (utr3_s, bar_y0), xmax - utr3_s, bh,
        color=COL_UTR, lw=0, zorder=6,
    ))

    # Boundary tick lines (full height)
    for xb in (cds_s, utr3_s):
        ax.axvline(xb, color=COL_BOUND, lw=0.9, ls='--', alpha=0.55, zorder=5)

    # Labels inside bar
    mid_utr5 = cds_s / 2
    mid_cds  = (cds_s + utr3_s) / 2
    mid_utr3 = (utr3_s + xmax) / 2
    bar_mid  = (bar_y0 + bar_y1) / 2

    fk = dict(fontsize=7.5, va='center', ha='center',
              fontweight='bold', color='white', zorder=7)

    if cds_s > 30:
        ax.text(mid_utr5, bar_mid, "5'UTR", **fk)
    ax.text(mid_cds, bar_mid, 'CDS', **fk)
    if xmax - utr3_s > 60:
        ax.text(mid_utr3, bar_mid, "3'UTR", **fk)

    # Boundary position ticks below x-axis
    for xb, lbl in ((cds_s, f'{cds_s}'), (utr3_s, f'{utr3_s}')):
        ax.text(xb, -0.07, lbl, ha='center', va='top',
                fontsize=6.5, color=COL_BOUND, zorder=5)


def main():
    df = pd.read_csv(HERE / 'shared_positions.tsv', sep='\t')

    fig = plt.figure(figsize=(10.5, 8.5))
    gs  = gridspec.GridSpec(3, 1, hspace=0.52)

    for i, gene in enumerate(GENES):
        ax  = fig.add_subplot(gs[i])
        sub = df[df['functional_gene'] == gene].sort_values('position')
        pos   = sub['position'].values
        flags = sub['shared_90pct'].values.astype(float)
        sm    = smooth(flags, SMOOTH)

        b    = BOUNDARIES[gene]
        xmax = pos[-1]

        # ── filled area ───────────────────────────────────────────────────────
        ax.fill_between(pos, 0, 1,    color=COL_UNIQUE, alpha=1.0, lw=0, zorder=1)
        ax.fill_between(pos, 0, sm,   color=COL_SHARED, alpha=0.85, lw=0, zorder=2)
        ax.plot(pos, sm, color=COL_LINE, lw=1.0, zorder=3)
        ax.axhline(0.5, color=COL_50, lw=0.8, ls='--', zorder=4, alpha=0.7)

        # ── annotation bar ────────────────────────────────────────────────────
        draw_annotation(ax, b, xmax)

        # ── overall percentage ────────────────────────────────────────────────
        pct = 100 * flags.mean()
        ax.text(0.985, 0.70,
                f'{pct:.0f}% non-unique',
                transform=ax.transAxes, ha='right', va='top',
                fontsize=9, color=COL_LINE, fontweight='bold')

        # ── title ─────────────────────────────────────────────────────────────
        ax.set_title(
            f'{gene}   ({PSEUDOGENES[gene]})',
            loc='left', fontsize=10.5, fontweight='bold', pad=14,
        )

        # ── axes formatting ───────────────────────────────────────────────────
        ax.set_xlim(0, xmax)
        ax.set_ylim(-0.04, 1.16)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'], fontsize=8)
        ax.set_ylabel('shared ≥90%\nwith pseudogene', fontsize=8, labelpad=4)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#aab7b8')
        ax.spines['bottom'].set_color('#aab7b8')
        ax.tick_params(axis='both', colors='#566573', labelsize=8)

        if gene != GENES[-1]:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel(
                'Position in MANE Select transcript (bp, 100 bp window start)',
                fontsize=9,
            )

    # ── legend ────────────────────────────────────────────────────────────────
    handles = [
        mpatches.Patch(color=COL_SHARED, alpha=0.85,
                       label='≥90% identical to ≥1 pseudogene  (100 bp window)'),
        mpatches.Patch(color=COL_UNIQUE,
                       label='unique region'),
        mpatches.Patch(color=COL_CDS,   label='CDS'),
        mpatches.Patch(color=COL_UTR,   label="5'/3' UTR"),
        plt.Line2D([0], [0], color=COL_50, lw=0.9, ls='--', label='50%'),
    ]
    fig.legend(handles=handles, loc='lower center', ncol=5,
               fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.03))

    fig.text(
        0.5, -0.07,
        'Smoothed with 80 bp rolling mean. GENCODE v49, MAFFT L-INS-i. '
        'Fractions are lower bounds — pseudogene transcript models are partial.',
        ha='center', fontsize=7, color='#717d7e',
    )

    for fmt in ('pdf', 'png'):
        out = HERE / f'hla_shared_fraction.{fmt}'
        fig.savefig(out, bbox_inches='tight', dpi=200)
        print(f'Saved → {out.name}')

    plt.show()


if __name__ == '__main__':
    main()
