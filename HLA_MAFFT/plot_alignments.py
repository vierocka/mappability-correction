"""
Visualise the three HLA MAFFT alignments.

Two output figures:
  Figure A — MANE mRNA  +  Genomic gene body
  Figure B — Deduplicated CDS (no cell-value annotations; wide cells for readability)

Each panel:
  - Left  : pairwise % identity heatmap
  - Right : per-position conservation track (sliding window = 100 bp)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from pathlib import Path

OUT    = Path(__file__).resolve().parent
WINDOW = 100   # sliding window for conservation track (bp) — matches read length

DEDUP_NOTE = (
    "\n\n\n\n\n"
    "Deduplicated: CDS sequences with identical nucleotide sequence were removed "
    "(seqkit rmdup -s). Multiple alleles and alternative loci are retained."
)

# ── Helpers ────────────────────────────────────────────────────────────────────
def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: ''.join(v).upper() for k, v in seqs.items()}

def count_seqs(path):
    if not path.exists():
        return 6
    return sum(1 for l in open(path) if l.startswith('>'))

def pairwise_identity(s1, s2):
    match = total = 0
    for a, b in zip(s1, s2):
        if a == '-' and b == '-': continue
        total += 1
        if a == b: match += 1
    return 100 * match / total if total else 0.0

def conservation_track(seqs_list):
    aln   = np.array([[c for c in s] for s in seqs_list])
    n_pos = aln.shape[1]
    cons  = np.zeros(n_pos)
    for i in range(n_pos):
        col     = aln[:, i]
        non_gap = col[col != '-']
        if len(non_gap) < 2: continue
        cons[i] = 1.0 if len(set(non_gap)) == 1 else 0.0
    kernel = np.ones(WINDOW) / WINDOW
    return np.convolve(cons, kernel, mode='same')

def draw_figure(alignments,
                annotate_cells=True,
                note=None,
                fig_width=20,
                fig_height_per_seq=0.35,
                col_width_ratios=None,
                font_scale=1.0):
    """
    alignments        : list of (title, path)
    annotate_cells    : print % values inside heatmap cells
    note              : (title_match, text) dimgray footnote for one panel
    fig_width         : total figure width in inches
    fig_height_per_seq: inches per sequence row (controls panel height)
    col_width_ratios  : [heatmap_weight, conservation_weight]
    font_scale        : multiplier applied to tick and annotation font sizes
    """
    if col_width_ratios is None:
        col_width_ratios = [1, 2.5]

    n_panels      = len(alignments)
    seq_counts    = [count_seqs(p) for _, p in alignments]
    height_ratios = [max(6.0, n * fig_height_per_seq) for n in seq_counts]
    total_height  = sum(height_ratios) + 1.5

    fig   = plt.figure(figsize=(fig_width, total_height))
    outer = gridspec.GridSpec(
        n_panels, 1, figure=fig,
        hspace=0.35,
        height_ratios=height_ratios)

    for panel_i, (title, aln_path) in enumerate(alignments):
        if not aln_path.exists():
            print(f"  SKIP (not found): {aln_path.name}")
            continue

        seqs  = read_fasta(aln_path)
        names = list(seqs.keys())
        n     = len(names)
        print(f"{title}: {n} sequences, alignment length {len(next(iter(seqs.values())))}")

        inner = gridspec.GridSpecFromSubplotSpec(
            1, 2, subplot_spec=outer[panel_i],
            width_ratios=col_width_ratios, wspace=0.1)

        # ── Left: pairwise identity heatmap ──────────────────────────────────
        ax_heat = fig.add_subplot(inner[0])
        mat = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                mat[i, j] = pairwise_identity(seqs[names[i]], seqs[names[j]])

        fs_tick = max(6.0, min(9.0, 60 / n)) * font_scale
        fs_ann  = max(5.5, min(7.5, 50 / n)) * font_scale

        im = ax_heat.imshow(mat, vmin=40, vmax=100, cmap='YlOrRd', aspect='auto')
        ax_heat.set_xticks(range(n))
        ax_heat.set_yticks(range(n))
        ax_heat.set_xticklabels(names, rotation=90, fontsize=fs_tick)
        ax_heat.set_yticklabels(names, fontsize=fs_tick)

        if annotate_cells:
            for i in range(n):
                for j in range(n):
                    ax_heat.text(j, i, f"{mat[i,j]:.0f}",
                                 ha='center', va='center', fontsize=fs_ann,
                                 color='black' if mat[i,j] < 85 else 'white')

        plt.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04, label='% identity')
        ax_heat.set_title(f"{title}\nPairwise % identity", fontsize=10, pad=8)

        if note and note[0] in title:
            ax_heat.text(
                0.5, -0.06, note[1],
                transform=ax_heat.transAxes,
                fontsize=7, color='dimgray',
                ha='center', va='top', style='italic')

        # ── Right: per-position conservation ─────────────────────────────────
        ax_cons = fig.add_subplot(inner[1])
        seq_list = list(seqs.values())
        aln_len  = len(seq_list[0])
        track    = conservation_track(seq_list)
        x        = np.arange(aln_len)

        ax_cons.fill_between(x, track, alpha=0.65, color='#4393c3', linewidth=0)
        ax_cons.set_xlim(0, aln_len)
        ax_cons.set_ylim(0, 1.05)
        ax_cons.set_xlabel("Alignment position (bp)", fontsize=9)
        ax_cons.set_ylabel("Conservation (sliding window)", fontsize=9)
        ax_cons.set_title(
            f"{title}\nPer-position conservation (window = {WINDOW} bp)",
            fontsize=10, pad=8)
        ax_cons.axhline(1.0, color='#d6604d', lw=0.7, linestyle='--',
                        label='fully conserved')
        ax_cons.legend(fontsize=8, loc='lower right')
        ax_cons.spines[['top', 'right']].set_visible(False)

    return fig

# ── Figure A: MANE mRNA + Genomic gene body ────────────────────────────────────
fig_a = draw_figure([
    ("MANE mRNA",
     OUT / "HLA_MANE.mafft.fna"),
    ("Genomic gene body (only canonical chromosomes)",
     OUT / "HLA_primaryChr.mafft.fna"),
], annotate_cells=True)

base_a = OUT / "HLA_alignments_A"
fig_a.savefig(str(base_a) + ".pdf", dpi=150, bbox_inches='tight')
fig_a.savefig(str(base_a) + ".png", dpi=150, bbox_inches='tight')
print(f"Saved: {base_a}.pdf / .png")
plt.close(fig_a)

# ── Figure B: Deduplicated CDS ────────────────────────────────────────────────
# 78 sequences — cells tripled in width relative to height:
#   fig_height_per_seq=0.18 → axes_h ≈ 14 in
#   fig_width=60            → heatmap axes ≈ 42 in wide → cells ~3× wider than tall
n_b      = count_seqs(OUT / "HLA_dedup_cds.mafft.fna")
fig_h_b  = n_b * 0.18 + 5

fig_b = draw_figure([
    ("Deduplicated CDS (also alternative)",
     OUT / "HLA_dedup_cds.mafft.fna"),
], annotate_cells=False,
   note=("Deduplicated", DEDUP_NOTE),
   fig_width=60,
   fig_height_per_seq=0.18,
   col_width_ratios=[3, 1],
   font_scale=2.0)

base_b = OUT / "HLA_alignments_B"
fig_b.savefig(str(base_b) + ".pdf", dpi=150, bbox_inches='tight', pad_inches=0.5)
fig_b.savefig(str(base_b) + ".png", dpi=150, bbox_inches='tight', pad_inches=0.5)
print(f"Saved: {base_b}.pdf / .png")
plt.close(fig_b)
