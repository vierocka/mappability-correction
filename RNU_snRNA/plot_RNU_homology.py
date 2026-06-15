#!/usr/bin/env python3
"""
Three-panel figure focused on RNU1 and RNU6 — the two snRNA families with UF=0.

Panel A: RNU1 family pairwise % identity (4×4 heatmap)
Panel B: RNU6 family pairwise % identity (6×6 heatmap, ATAC variant highlighted)
Panel C: Uniqueness factor bar chart — RNU1-1/2/3/4 and RNU6-1/2/7/8/9
         across four MANE-reference simulation settings
"""

import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

HERE    = Path(__file__).resolve().parent
RESULTS = HERE.parent / "mappability" / "results"

# ── load per-family MAFFT alignment ──────────────────────────────────────────
def load_mafft(path):
    seqs, names = {}, []
    cur = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                cur = line[1:].split('|')[0]
                names.append(cur)
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line)
    return names, {k: ''.join(v) for k, v in seqs.items()}

def pct_identity(a, b):
    pairs = [(x, y) for x, y in zip(a, b) if x != '-' and y != '-']
    if not pairs:
        return 0.0
    return 100 * sum(x == y for x, y in pairs) / len(pairs)

def identity_matrix(names, aligned):
    n = len(names)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            mat[i, j] = pct_identity(aligned[names[i]], aligned[names[j]])
    return mat

rnu1_names, rnu1_aln = load_mafft(HERE / 'RNU1_family.mafft.fna')
rnu6_names, rnu6_aln = load_mafft(HERE / 'RNU6_family.mafft.fna')

rnu1_mat = identity_matrix(rnu1_names, rnu1_aln)
rnu6_mat = identity_matrix(rnu6_names, rnu6_aln)

# ── UF data for panel C ───────────────────────────────────────────────────────
# IDs from RNU_MANE.fna
focus_genes = {
    'ENST00000383925.1': 'RNU1-1',
    'ENST00000384278.1': 'RNU1-2',
    'ENST00000384782.1': 'RNU1-3',
    'ENST00000384659.1': 'RNU1-4',
    'ENST00000383898.1': 'RNU6-1',
    'ENST00000384627.1': 'RNU6-2',
    'ENST00000364784.1': 'RNU6-7',
    'ENST00000365467.1': 'RNU6-8',
    'ENST00000384776.1': 'RNU6-9',
}

settings = {
    'SE 75 bp':            'transcript_uniqueness_factors_MANE_RNA_SE_L75bp.tsv',
    'SE 100 bp':           'transcript_uniqueness_factors_MANE_RNA_SE_L100bp.tsv',
    'exon-flank PE 100 bp': 'transcript_uniqueness_factors_MANE_exonflank_PE_L100bp.tsv',
}
# SE 200 bp is excluded from Panel C: all 9 genes (107–164 bp) are shorter than
# 200 bp and cannot generate reads, producing only × markers with no signal.

uf_data = {gene: {} for gene in focus_genes.values()}
for label, fname in settings.items():
    path = RESULTS / fname
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row['transcript_id'] in focus_genes:
                gene = focus_genes[row['transcript_id']]
                uf_data[gene][label] = float(row['uniqueness_factor'])
    for gene in uf_data:
        if label not in uf_data[gene]:
            uf_data[gene][label] = np.nan   # too short

# ── figure layout ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 5))
gs  = fig.add_gridspec(1, 3, width_ratios=[0.9, 1.35, 1.6], wspace=0.38)
ax_rnu1 = fig.add_subplot(gs[0])
ax_rnu6 = fig.add_subplot(gs[1])
ax_bar  = fig.add_subplot(gs[2])

# ── helper: draw heatmap ──────────────────────────────────────────────────────
def draw_heatmap(ax, names, mat, title, vmin=90, vmax=100, annotate_thresh=90):
    n   = len(names)
    # colour palette: white (identical off-diagonal looks like white if 100%)
    import matplotlib.colors as mcolors
    cmap = plt.cm.YlOrRd
    im   = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(names, fontsize=8)
    for i in range(n):
        for j in range(n):
            val = mat[i, j]
            if i != j and val >= annotate_thresh:
                ax.text(j, i, f'{val:.1f}',
                        ha='center', va='center', fontsize=7,
                        color='white' if val >= 99 else '#333333',
                        fontweight='bold')
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label('Pairwise % identity', fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    ax.set_title(title, fontsize=9, fontweight='bold', loc='left', pad=7)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    return im

draw_heatmap(ax_rnu1, rnu1_names, rnu1_mat,
             'A   RNU1 family\n    (MANE Select, 4 loci)',
             vmin=98, vmax=100)

draw_heatmap(ax_rnu6, rnu6_names, rnu6_mat,
             'B   RNU6 family\n    (MANE Select, 6 loci)',
             vmin=70, vmax=100, annotate_thresh=80)

# ── Panel C: UF heatmap (genes × settings) ───────────────────────────────────
# A bar chart is misleading when all values are 0 (invisible bars).
# A colour-annotated heatmap makes the zero-signal message explicit.
gene_order = ['RNU1-1','RNU1-2','RNU1-3','RNU1-4',
              'RNU6-1','RNU6-2','RNU6-7','RNU6-8','RNU6-9']
set_labels = list(settings.keys())   # SE 75, SE 100, exflank PE 100
n_genes    = len(gene_order)
n_sets     = len(set_labels)

# Build matrix: rows = genes, cols = settings
mat_c = np.zeros((n_genes, n_sets))
for gi, gene in enumerate(gene_order):
    for si, label in enumerate(set_labels):
        v = uf_data[gene].get(label, np.nan)
        mat_c[gi, si] = 0.0 if np.isnan(v) else v

# Use a green-white-red diverging palette: 0=red, 1=green
cmap_c = plt.cm.RdYlGn
im_c   = ax_bar.imshow(mat_c, cmap=cmap_c, vmin=0, vmax=1, aspect='auto')

# Annotate every cell with the numeric value
for gi in range(n_genes):
    for si in range(n_sets):
        val = mat_c[gi, si]
        ax_bar.text(si, gi, f'{val:.3f}', ha='center', va='center',
                    fontsize=8, fontweight='bold',
                    color='white' if val < 0.25 or val > 0.75 else '#333333')

ax_bar.set_xticks(range(n_sets))
ax_bar.set_xticklabels(set_labels, fontsize=8)
ax_bar.set_yticks(range(n_genes))
ax_bar.set_yticklabels(gene_order, fontsize=8)

# Horizontal separator between RNU1 and RNU6 groups
ax_bar.axhline(3.5, color='white', linewidth=2.5)

# Family labels on the right
ax_bar.annotate('U1 snRNA\n(4 identical copies)',
                xy=(1.03, 1 - (1.5 / n_genes)),
                xycoords='axes fraction', fontsize=7.5, va='center',
                color='#2166ac', style='italic')
ax_bar.annotate('U6 snRNA\n(5 near-identical copies)',
                xy=(1.03, 1 - (6.5 / n_genes)),
                xycoords='axes fraction', fontsize=7.5, va='center',
                color='#d6604d', style='italic')

cb_c = fig.colorbar(im_c, ax=ax_bar, orientation='horizontal',
                    location='bottom', fraction=0.04, pad=0.22, shrink=0.85)
cb_c.set_label('Uniqueness factor', fontsize=7.5, labelpad=3)
cb_c.ax.tick_params(labelsize=7)

ax_bar.set_title('C   Uniqueness factor — 3 simulation settings\n'
                 '    (MANE reference; reads from MANE RNA)',
                 fontsize=9, fontweight='bold', loc='left', pad=7)
for spine in ax_bar.spines.values():
    spine.set_visible(False)

fig.suptitle('RNU1 and RNU6 snRNAs: high sequence identity causes zero uniqueness factor in MANE reference',
             fontsize=10.5, fontweight='bold', y=1.02)

for ext in ('pdf', 'png'):
    out = HERE / f'RNU_homology.{ext}'
    fig.savefig(out, dpi=180, bbox_inches='tight')
    print(f'Saved → {out}')

plt.close()
