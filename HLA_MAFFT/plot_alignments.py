"""
Visualise the three HLA MAFFT alignments.

For each alignment:
  - Pairwise % identity heatmap
  - Per-position conservation track (fraction of non-gap identical columns)

Output: HLA_alignments.pdf / .png
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from pathlib import Path

OUT = Path("/home/veve/Dropbox/Self-Nonself/Reference/mappability-correction/HLA_MAFFT")

ALIGNMENTS = [
    ("MANE mRNA",         OUT / "HLA_MANE.mafft.fna"),
    ("Dedup CDS",         OUT / "HLA_dedup_cds.mafft.fna"),
    ("Primary chr gene body (no alt contigs)", OUT / "HLA_primaryChr.mafft.fna"),
]

# ── Parsers ────────────────────────────────────────────────────────────────────
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

def pairwise_identity(s1, s2):
    """% identity over aligned (non-double-gap) columns."""
    match = total = 0
    for a, b in zip(s1, s2):
        if a == '-' and b == '-': continue
        total += 1
        if a == b: match += 1
    return 100 * match / total if total else 0.0

def conservation_track(seqs_list, window=50):
    """Fraction of positions in a sliding window that are fully conserved."""
    aln = np.array([[c for c in s] for s in seqs_list])
    n_pos = aln.shape[1]
    cons  = np.zeros(n_pos)
    for i in range(n_pos):
        col = aln[:, i]
        non_gap = col[col != '-']
        if len(non_gap) < 2: continue
        cons[i] = 1.0 if len(set(non_gap)) == 1 else 0.0
    # Sliding window average
    kernel = np.ones(window) / window
    return np.convolve(cons, kernel, mode='same')

# ── Figure ─────────────────────────────────────────────────────────────────────
n_panels = len(ALIGNMENTS)
fig = plt.figure(figsize=(18, 6 * n_panels))
outer = gridspec.GridSpec(n_panels, 1, figure=fig, hspace=0.55)

for panel_i, (title, aln_path) in enumerate(ALIGNMENTS):
    if not aln_path.exists():
        print(f"  SKIP (not found): {aln_path.name}")
        continue

    seqs = read_fasta(aln_path)
    names = list(seqs.keys())
    n = len(names)
    print(f"{title}: {n} sequences, alignment length {len(next(iter(seqs.values())))}")

    inner = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=outer[panel_i],
        width_ratios=[1, 2.5], wspace=0.35)

    # ── Left: pairwise identity heatmap ──────────────────────────────────────
    ax_heat = fig.add_subplot(inner[0])
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            mat[i, j] = pairwise_identity(seqs[names[i]], seqs[names[j]])

    im = ax_heat.imshow(mat, vmin=40, vmax=100, cmap='YlOrRd', aspect='auto')
    ax_heat.set_xticks(range(n))
    ax_heat.set_yticks(range(n))
    ax_heat.set_xticklabels(names, rotation=90, fontsize=7)
    ax_heat.set_yticklabels(names, fontsize=7)
    for i in range(n):
        for j in range(n):
            ax_heat.text(j, i, f"{mat[i,j]:.0f}",
                         ha='center', va='center', fontsize=5.5,
                         color='black' if mat[i,j] < 85 else 'white')
    plt.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04, label='% identity')
    ax_heat.set_title(f"{title}\nPairwise % identity", fontsize=10, pad=8)

    # ── Right: per-position conservation ─────────────────────────────────────
    ax_cons = fig.add_subplot(inner[1])
    seq_list = list(seqs.values())
    aln_len  = len(seq_list[0])
    track    = conservation_track(seq_list, window=max(1, aln_len // 200))
    x        = np.arange(aln_len)

    ax_cons.fill_between(x, track, alpha=0.65, color='#4393c3', linewidth=0)
    ax_cons.set_xlim(0, aln_len)
    ax_cons.set_ylim(0, 1.05)
    ax_cons.set_xlabel("Alignment position (bp)", fontsize=9)
    ax_cons.set_ylabel("Conservation (sliding window)", fontsize=9)
    ax_cons.set_title(f"{title}\nPer-position conservation (window={max(1, aln_len//200)} bp)",
                      fontsize=10, pad=8)
    ax_cons.axhline(1.0, color='#d6604d', lw=0.7, linestyle='--', label='fully conserved')
    ax_cons.legend(fontsize=8, loc='lower right')
    ax_cons.spines[['top', 'right']].set_visible(False)

base = OUT / "HLA_alignments"
fig.savefig(str(base) + ".pdf", dpi=150, bbox_inches='tight')
fig.savefig(str(base) + ".png", dpi=150, bbox_inches='tight')
print(f"\nSaved: {base}.pdf / .png")
