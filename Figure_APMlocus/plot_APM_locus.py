"""
APM locus — MHC class I region, GRCh38.p14
Two-panel figure:
  Panel A : TAP/PSM/TAPBP cluster and flanking genes  (~32.75 – 33.40 Mb)
            All 50 gene/pseudogene features from NCBI GFF included.
  Panel B : HLA class I cluster    (~29.65 – 31.60 Mb)

Multi-track layout: genes on the same strand that would visually overlap
are placed on separate horizontal tracks.

Dotted yellow vertical lines mark the genomic position where any two
annotated features overlap (regardless of strand).

Source: GCF_000001405.40_GRCh38.p14_genomic.gff  (NCBI GRCh38.p14)
Snapshot: APM_locus_GFF_snapshot.tsv
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np

# ── Colour palette ─────────────────────────────────────────────────────────────
C_POS     = '#2166ac'   # blue  – protein-coding, + strand
C_NEG     = '#d6604d'   # red   – protein-coding, - strand
C_PSEUDO  = '#999999'   # grey  – pseudogene (outline only)
C_AXIS    = '#333333'
C_OVERLAP = '#BA55D3'   # medium orchid – overlap marker

GENE_H    = 0.26   # box height
TRACK_SEP = 0.60   # vertical distance between tracks

# ── Panel A: all 50 features in 32.75–33.40 Mb from NCBI GFF ──────────────────
# (start, end, strand, label, is_pseudo)
PANEL_A_GENES = [
    (32756098, 32763532, '-', 'HLA-DQB2',     False),
    (32812763, 32817002, '-', 'HLA-DOB',      False),
    (32821831, 32838739, '-', 'TAP2',         False),
    (32840717, 32844679, '-', 'PSMB8',        False),
    (32844086, 32846500, '+', 'PSMB8-AS1',    False),
    (32845209, 32853704, '-', 'TAP1',         False),
    (32854192, 32859851, '+', 'PSMB9',        False),
    (32876478, 32880074, '-', 'PPP1R2P1',     True),
    (32894176, 32903758, '+', 'LOC100294145', False),
    (32896402, 32896489, '+', 'HLA-Z',        True),
    (32934636, 32941028, '-', 'HLA-DMB',      False),
    (32948618, 32953097, '-', 'HLA-DMA',      False),
    (32963293, 32968521, '-', 'LOC124901302', False),
    (32968594, 32981505, '+', 'BRD2',         False),
    (33004182, 33009591, '-', 'HLA-DOA',      False),
    (33064569, 33080748, '-', 'HLA-DPA1',     False),
    (33075990, 33089696, '+', 'HLA-DPB1',     False),
    (33079398, 33079909, '+', 'RPL32P1',      True),
    (33091482, 33093314, '-', 'HLA-DPA2',     True),
    (33103693, 33107144, '-', 'COL11A2P1',    True),
    (33112516, 33129113, '+', 'HLA-DPB2',     True),
    (33130747, 33143434, '-', 'LOC105375021', False),
    (33131197, 33131343, '-', 'HLA-DPA3',     True),
    (33143324, 33148815, '+', 'HCG24',        False),
    (33162694, 33193519, '-', 'COL11A2',      False),
    (33193588, 33200853, '-', 'RXRB',         False),
    (33199601, 33199692, '+', 'RNY4P10',      True),
    (33200867, 33204437, '+', 'SLC39A7',      False),
    (33204655, 33206831, '+', 'HSD17B8',      False),
    (33207835, 33207944, '+', 'MIR219A1',     False),
    (33208500, 33212716, '+', 'RING1',        False),
    (33215574, 33216183, '+', 'ZNF70P1',      True),
    (33246075, 33246873, '-', 'LOC105375022', False),
    (33249536, 33254890, '+', 'HCG25',        False),
    (33250272, 33271965, '-', 'VPS52',        False),
    (33272075, 33276511, '+', 'RPS18',        False),
    (33277123, 33278825, '+', 'B3GALT4',      False),
    (33279108, 33289239, '-', 'WDR46',        False),
    (33287227, 33287289, '-', 'MIR6873',      False),
    (33289197, 33290934, '+', 'PFDN6',        False),
    (33290245, 33290325, '+', 'MIR6834',      False),
    (33291654, 33301612, '-', 'RGL2',         False),
    (33299694, 33314078, '-', 'TAPBP',        False),
    (33314418, 33317942, '-', 'ZBTB22',       False),
    (33318558, 33322959, '-', 'DAXX',         False),
    (33322674, 33323647, '+', 'LOC124901303', False),
    (33323628, 33329279, '-', 'SMIM40',       False),
    (33338963, 33339711, '+', 'MYL12BP3',     True),
    (33364737, 33366362, '-', 'LYPLA2P1',     True),
    (33389356, 33389478, '-', 'RPL35AP4',     True),
]

# ── Panel B: HLA class I cluster ───────────────────────────────────────────────
PANEL_B_GENES = [
    # Protein-coding
    (29723434, 29738532, '+', 'HLA-F', False),
    (29826474, 29831021, '+', 'HLA-G', False),
    (29942532, 29945870, '+', 'HLA-A', False),
    (30489509, 30494194, '+', 'HLA-E', False),
    (31268749, 31272092, '-', 'HLA-C', False),
    (31353875, 31357179, '-', 'HLA-B', False),
    (31400711, 31415315, '+', 'MICA',  False),
    (31494918, 31511124, '+', 'MICB',  False),
    # Pseudogenes – HLA class I
    (29791906, 29797807, '+', 'HLA-V', True),
    (29800044, 29803079, '+', 'HLA-P', True),
    (29887573, 29891079, '+', 'HLA-H', True),
    (29896443, 29898947, '+', 'HLA-T', True),
    (29926659, 29929825, '+', 'HLA-K', True),
    (29933764, 29934880, '+', 'HLA-U', True),
    (29955834, 29959058, '+', 'HLA-W', True),
    (30005971, 30009956, '+', 'HLA-J', True),
    (30259562, 30266951, '+', 'HLA-L', True),
    (30351074, 30352038, '+', 'HLA-N', True),
    # Pseudogenes – MIC
    (29741556, 29745784, '-', 'MICE',  True),
    (29812390, 29812692, '-', 'MICG',  True),
    (29852187, 29854052, '-', 'MICF',  True),
    (29970372, 29975745, '-', 'MICD',  True),
]

# ── Track assignment ───────────────────────────────────────────────────────────
def assign_tracks(genes, gap_frac=0.008, region_span=1):
    gap = region_span * gap_frac
    pos_tracks, neg_tracks = [], []
    result = {}
    for (start, end, strand, label, _) in sorted(genes, key=lambda g: g[0]):
        pool = pos_tracks if strand == '+' else neg_tracks
        placed = False
        for i, t_end in enumerate(pool):
            if start >= t_end + gap:
                pool[i] = end
                result[label] = i
                placed = True
                break
        if not placed:
            pool.append(end)
            result[label] = len(pool) - 1
    return result

# ── Overlap detection ──────────────────────────────────────────────────────────
def find_overlaps(genes):
    """Return list of (x_overlap_start, label1, label2) for all overlapping pairs."""
    overlaps = []
    by_start = sorted(genes, key=lambda g: g[0])
    for i in range(len(by_start)):
        s1, e1, _, l1, _ = by_start[i]
        for j in range(i + 1, len(by_start)):
            s2, e2, _, l2, _ = by_start[j]
            if s2 >= e1:
                break
            overlaps.append((s2, l1, l2))
    return overlaps

# ── Drawing helpers ────────────────────────────────────────────────────────────
def gene_color(strand, is_pseudo):
    if is_pseudo: return C_PSEUDO
    return C_POS if strand == '+' else C_NEG

def draw_gene(ax, start, end, strand, label, is_pseudo, track, region_span,
              label_fontsize=7.5, min_disp_frac=0.004):
    color  = gene_color(strand, is_pseudo)
    filled = not is_pseudo
    y_base = (track + 1) * TRACK_SEP * (1 if strand == '+' else -1)

    # Ensure tiny genes are visible
    min_w  = region_span * min_disp_frac
    raw_w  = end - start
    disp_w = max(raw_w, min_w)
    cx     = (start + end) / 2
    ds     = cx - disp_w / 2   # display start

    arrow_w = min(region_span * 0.015, disp_w * 0.40)
    body_w  = max(disp_w - arrow_w, disp_w * 0.02)

    fc = color if filled else 'white'
    if strand == '+':
        ax.add_patch(plt.Rectangle((ds, y_base - GENE_H/2), body_w, GENE_H,
                                   fc=fc, ec=color, lw=0.9, zorder=3))
        ax.fill([ds + body_w, ds + body_w, ds + disp_w],
                [y_base - GENE_H/2, y_base + GENE_H/2, y_base],
                fc=fc, ec=color, lw=0.9, zorder=3)
    else:
        ax.add_patch(plt.Rectangle((ds + arrow_w, y_base - GENE_H/2), body_w, GENE_H,
                                   fc=fc, ec=color, lw=0.9, zorder=3))
        ax.fill([ds + arrow_w, ds + arrow_w, ds],
                [y_base - GENE_H/2, y_base + GENE_H/2, y_base],
                fc=fc, ec=color, lw=0.9, zorder=3)

    # Vertical stem + rotated label
    x_mid    = cx
    stem_gap = 0.05
    if strand == '+':
        stem_y0 = y_base + GENE_H / 2
        stem_y1 = stem_y0 + stem_gap + track * 0.14
        va, ha  = 'bottom', 'center'
    else:
        stem_y0 = y_base - GENE_H / 2
        stem_y1 = stem_y0 - stem_gap - track * 0.14
        va, ha  = 'top', 'center'

    ax.plot([x_mid, x_mid], [stem_y0, stem_y1],
            color='#bbbbbb', lw=0.5, zorder=2)
    style = 'italic' if is_pseudo else 'normal'
    ax.text(x_mid, stem_y1, label,
            ha=ha, va=va, fontsize=label_fontsize,
            fontstyle=style, rotation=90, rotation_mode='anchor',
            color='#555555' if is_pseudo else '#111111', clip_on=True)


def draw_panel(ax, genes, region_start, region_end, panel_label,
               label_fontsize=7.5, gap_frac=0.008, min_disp_frac=0.004):
    span   = region_end - region_start
    tracks = assign_tracks(genes, gap_frac=gap_frac, region_span=span)

    n_pos = max((tracks[l] for (_, _, s, l, _) in genes if s == '+'), default=-1) + 1
    n_neg = max((tracks[l] for (_, _, s, l, _) in genes if s == '-'), default=-1) + 1

    # Draw gene bodies and labels
    for g in genes:
        draw_gene(ax, *g[:4], g[4], tracks[g[3]], span,
                  label_fontsize, min_disp_frac)

    # Overlap markers — dotted yellow vertical lines
    track_y = {l: (tracks[l] + 1) * TRACK_SEP * (1 if s == '+' else -1)
                for (_, _, s, l, _) in genes}
    for (x, l1, l2) in find_overlaps(genes):
        y1   = track_y[l1]
        y2   = track_y[l2]
        y_lo = min(y1, y2) - GENE_H / 2
        y_hi = max(y1, y2) + GENE_H / 2
        ax.vlines(x, y_lo, y_hi,
                  colors=C_OVERLAP, linewidths=1.5,
                  linestyles='dotted', zorder=6)

    # Chromosome backbone
    ax.hlines(0, region_start, region_end, colors=C_AXIS, linewidths=1.8, zorder=2)

    # Strand labels
    ax.text(region_start - span * 0.013, TRACK_SEP * 0.5, "5'→3'",
            ha='right', va='center', fontsize=8.5, color=C_POS, fontweight='bold')
    ax.text(region_start - span * 0.013, -TRACK_SEP * 0.5, "3'→5'",
            ha='right', va='center', fontsize=8.5, color=C_NEG, fontweight='bold')

    # Scale bar
    sb_len = 200_000 if span < 800_000 else 500_000
    sb_x   = region_end - sb_len * 1.3
    y_sb   = -(n_neg + 1.6) * TRACK_SEP
    ax.hlines(y_sb, sb_x, sb_x + sb_len, colors=C_AXIS, linewidths=2)
    ax.vlines([sb_x, sb_x + sb_len], y_sb - 0.06, y_sb + 0.06,
              colors=C_AXIS, linewidths=1.5)
    ax.text(sb_x + sb_len / 2, y_sb - 0.13, f"{sb_len // 1000} kb",
            ha='center', va='top', fontsize=7.5)

    # Position ticks
    ticks = np.linspace(region_start, region_end, 5)
    for t in ticks:
        ax.vlines(t, -0.09, 0.09, colors='#999999', linewidths=0.8)
        ax.text(t, y_sb - 0.38, f"{t / 1e6:.2f}",
                ha='center', va='top', fontsize=6.5, color='#555555')
    ax.text((region_start + region_end) / 2, y_sb - 0.65,
            "chr6 position (Mb)",
            ha='center', va='top', fontsize=7.5, color='#444444')

    # Panel label
    ax.text(region_start - span * 0.013, (n_pos + 1.8) * TRACK_SEP,
            panel_label, fontsize=13, fontweight='bold', va='bottom')

    ax.set_xlim(region_start - span * 0.016, region_end + span * 0.005)
    ax.set_ylim(y_sb - 0.85, (n_pos + 2.8) * TRACK_SEP)
    ax.axis('off')


# ── Build figure ───────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(17, 10))
gs  = gridspec.GridSpec(2, 1, figure=fig, hspace=0.25,
                        height_ratios=[1.4, 1.0])

ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1])

draw_panel(ax_a, PANEL_A_GENES,
           region_start=32_750_000, region_end=33_400_000,
           panel_label='A', label_fontsize=6.5, gap_frac=0.003,
           min_disp_frac=0.004)

draw_panel(ax_b, PANEL_B_GENES,
           region_start=29_650_000, region_end=31_600_000,
           panel_label='B', label_fontsize=7.5, gap_frac=0.008,
           min_disp_frac=0.003)

# ── Legend ─────────────────────────────────────────────────────────────────────
legend_elements = [
    mpatches.Patch(fc=C_POS,     ec=C_POS,     label="Protein-coding, 5'→3' (+)"),
    mpatches.Patch(fc=C_NEG,     ec=C_NEG,     label="Protein-coding, 3'→5' (−)"),
    mpatches.Patch(fc='white',   ec=C_PSEUDO,  label='Pseudogene'),
    mpatches.Patch(fc=C_OVERLAP, ec=C_OVERLAP, label='Gene overlap'),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=4,
           fontsize=9, framealpha=0.9, bbox_to_anchor=(0.5, 0.005))

base = 'mappability-correction/Figure_APMlocus/APM_locus_chr6'
fig.savefig(base + '.pdf', dpi=200, bbox_inches='tight')
fig.savefig(base + '.png', dpi=200, bbox_inches='tight')
print(f"Saved: {base}.pdf")
print(f"Saved: {base}.png")
