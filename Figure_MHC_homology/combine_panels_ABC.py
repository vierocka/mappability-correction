#!/usr/bin/env python3
"""
Produces two figures:
  Figure_ABC_vertical.pdf/png  — A, B, C stacked (legend under B shrunk + spacer)
  Figure_ABC_vchrom.pdf/png    — A (vertical chr line) | B (vertical chr line) | C
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from pathlib import Path

FS    = 2.0   # panel C font scale
FS_AB = 3.0   # panels A and B font scale
FS_C  = FS * 1.25   # = 2.5, panel C 1.25x font boost

DIR_MAFFT = Path(__file__).parent / 'HLApseudo_MAFFT'

# ── Colours ────────────────────────────────────────────────────────────────────
C_POS     = '#2166ac'
C_NEG     = '#d6604d'
C_PSEUDO  = '#999999'
C_AXIS    = '#333333'
C_OVERLAP = '#BA55D3'

GENE_H    = 0.26
TRACK_SEP = 0.60

# ── Gene lists ─────────────────────────────────────────────────────────────────
PANEL_A_GENES = [
    (32756098,32763532,'-','HLA-DQB2',False),(32812763,32817002,'-','HLA-DOB',False),
    (32821831,32838739,'-','TAP2',False),(32840717,32844679,'-','PSMB8',False),
    (32844086,32846500,'+','PSMB8-AS1',False),(32845209,32853704,'-','TAP1',False),
    (32854192,32859851,'+','PSMB9',False),(32876478,32880074,'-','PPP1R2P1',True),
    (32894176,32903758,'+','LOC100294145',False),(32896402,32896489,'+','HLA-Z',True),
    (32934636,32941028,'-','HLA-DMB',False),(32948618,32953097,'-','HLA-DMA',False),
    (32963293,32968521,'-','LOC124901302',False),(32968594,32981505,'+','BRD2',False),
    (33004182,33009591,'-','HLA-DOA',False),(33064569,33080748,'-','HLA-DPA1',False),
    (33075990,33089696,'+','HLA-DPB1',False),(33079398,33079909,'+','RPL32P1',True),
    (33091482,33093314,'-','HLA-DPA2',True),(33103693,33107144,'-','COL11A2P1',True),
    (33112516,33129113,'+','HLA-DPB2',True),(33130747,33143434,'-','LOC105375021',False),
    (33131197,33131343,'-','HLA-DPA3',True),(33143324,33148815,'+','HCG24',False),
    (33162694,33193519,'-','COL11A2',False),(33193588,33200853,'-','RXRB',False),
    (33199601,33199692,'+','RNY4P10',True),(33200867,33204437,'+','SLC39A7',False),
    (33204655,33206831,'+','HSD17B8',False),(33207835,33207944,'+','MIR219A1',False),
    (33208500,33212716,'+','RING1',False),(33215574,33216183,'+','ZNF70P1',True),
    (33246075,33246873,'-','LOC105375022',False),(33249536,33254890,'+','HCG25',False),
    (33250272,33271965,'-','VPS52',False),(33272075,33276511,'+','RPS18',False),
    (33277123,33278825,'+','B3GALT4',False),(33279108,33289239,'-','WDR46',False),
    (33287227,33287289,'-','MIR6873',False),(33289197,33290934,'+','PFDN6',False),
    (33290245,33290325,'+','MIR6834',False),(33291654,33301612,'-','RGL2',False),
    (33299694,33314078,'-','TAPBP',False),(33314418,33317942,'-','ZBTB22',False),
    (33318558,33322959,'-','DAXX',False),(33322674,33323647,'+','LOC124901303',False),
    (33323628,33329279,'-','SMIM40',False),(33338963,33339711,'+','MYL12BP3',True),
    (33364737,33366362,'-','LYPLA2P1',True),(33389356,33389478,'-','RPL35AP4',True),
]

PANEL_B_GENES = [
    (29723434,29738532,'+','HLA-F',False),(29826474,29831021,'+','HLA-G',False),
    (29942532,29945870,'+','HLA-A',False),(30489509,30494194,'+','HLA-E',False),
    (31268749,31272092,'-','HLA-C',False),(31353875,31357179,'-','HLA-B',False),
    (31400711,31415315,'+','MICA',False),(31494918,31511124,'+','MICB',False),
    (29791906,29797807,'+','HLA-V',True),(29800044,29803079,'+','HLA-P',True),
    (29887573,29891079,'+','HLA-H',True),(29896443,29898947,'+','HLA-T',True),
    (29926659,29929825,'+','HLA-K',True),(29933764,29934880,'+','HLA-U',True),
    (29955834,29959058,'+','HLA-W',True),(30005971,30009956,'+','HLA-J',True),
    (30259562,30266951,'+','HLA-L',True),(30351074,30352038,'+','HLA-N',True),
    (29741556,29745784,'-','MICE',True),(29812390,29812692,'-','MICG',True),
    (29852187,29854052,'-','MICF',True),(29970372,29975745,'-','MICD',True),
]

# ── Helpers ────────────────────────────────────────────────────────────────────
def assign_tracks(genes, gap_frac=0.008, region_span=1):
    gap = region_span * gap_frac
    pos_tracks, neg_tracks = [], []
    result = {}
    for (s, e, strand, label, _) in sorted(genes, key=lambda g: g[0]):
        pool = pos_tracks if strand == '+' else neg_tracks
        placed = False
        for i, t_end in enumerate(pool):
            if s >= t_end + gap:
                pool[i] = e; result[label] = i; placed = True; break
        if not placed:
            pool.append(e); result[label] = len(pool) - 1
    return result

def find_overlaps(genes):
    overlaps = []
    by_start = sorted(genes, key=lambda g: g[0])
    for i in range(len(by_start)):
        s1,e1,_,l1,_ = by_start[i]
        for j in range(i+1, len(by_start)):
            s2,e2,_,l2,_ = by_start[j]
            if s2 >= e1: break
            overlaps.append((s2, l1, l2))
    return overlaps

def gene_color(strand, is_pseudo):
    if is_pseudo: return C_PSEUDO
    return C_POS if strand == '+' else C_NEG

def locus_legend_handles():
    return [
        mpatches.Patch(fc=C_POS,   ec=C_POS,   label="Protein-coding, 5'→3' (+)"),
        mpatches.Patch(fc=C_NEG,   ec=C_NEG,   label="Protein-coding, 3'→5' (−)"),
        mpatches.Patch(fc='white', ec=C_PSEUDO, label='Pseudogene'),
        mpatches.Patch(fc=C_OVERLAP,ec=C_OVERLAP,label='Gene overlap'),
    ]

# ── Horizontal locus panel (original orientation) ──────────────────────────────
def draw_gene_h(ax, start, end, strand, label, is_pseudo, track, span,
                label_fontsize, min_disp_frac=0.004):
    color = gene_color(strand, is_pseudo)
    fc    = color if not is_pseudo else 'white'
    yb    = (track + 1) * TRACK_SEP * (1 if strand == '+' else -1)
    dw    = max(end - start, span * min_disp_frac)
    cx    = (start + end) / 2
    ds    = cx - dw / 2
    aw    = min(span * 0.015, dw * 0.40)
    bw    = max(dw - aw, dw * 0.02)
    if strand == '+':
        ax.add_patch(plt.Rectangle((ds, yb-GENE_H/2), bw, GENE_H, fc=fc, ec=color, lw=0.9, zorder=3))
        ax.fill([ds+bw,ds+bw,ds+dw],[yb-GENE_H/2,yb+GENE_H/2,yb],fc=fc,ec=color,lw=0.9,zorder=3)
    else:
        ax.add_patch(plt.Rectangle((ds+aw, yb-GENE_H/2), bw, GENE_H, fc=fc, ec=color, lw=0.9, zorder=3))
        ax.fill([ds+aw,ds+aw,ds],[yb-GENE_H/2,yb+GENE_H/2,yb],fc=fc,ec=color,lw=0.9,zorder=3)
    sg = 0.05
    if strand == '+':
        sy0 = yb+GENE_H/2; sy1 = sy0+sg+track*0.14; va,ha = 'bottom','center'
    else:
        sy0 = yb-GENE_H/2; sy1 = sy0-sg-track*0.14; va,ha = 'top','center'
    ax.plot([cx,cx],[sy0,sy1],color='#bbbbbb',lw=0.5,zorder=2)
    ax.text(cx, sy1, label, ha=ha, va=va, fontsize=label_fontsize,
            fontstyle='italic' if is_pseudo else 'normal',
            rotation=90, rotation_mode='anchor',
            color='#555555' if is_pseudo else '#111111', clip_on=True)

def draw_panel_h(ax, genes, r0, r1, panel_label, label_fontsize, fs,
                 gap_frac=0.008, min_disp_frac=0.004):
    span   = r1 - r0
    tracks = assign_tracks(genes, gap_frac=gap_frac, region_span=span)
    n_pos  = max((tracks[l] for (_,_,s,l,_) in genes if s=='+'), default=-1)+1
    n_neg  = max((tracks[l] for (_,_,s,l,_) in genes if s=='-'), default=-1)+1
    for g in genes:
        draw_gene_h(ax, *g[:4], g[4], tracks[g[3]], span, label_fontsize, min_disp_frac)
    track_y = {l:(tracks[l]+1)*TRACK_SEP*(1 if s=='+' else -1) for (_,_,s,l,_) in genes}
    for (x,l1,l2) in find_overlaps(genes):
        y1,y2 = track_y[l1],track_y[l2]
        ax.vlines(x,min(y1,y2)-GENE_H/2,max(y1,y2)+GENE_H/2,
                  colors=C_OVERLAP,lw=1.5,linestyles='dotted',zorder=6)
    ax.hlines(0,r0,r1,colors=C_AXIS,lw=1.8,zorder=2)
    ax.text(r0-span*0.013, TRACK_SEP*0.5,"5'→3'",ha='right',va='center',fontsize=8.5*fs,color=C_POS,fontweight='bold')
    ax.text(r0-span*0.013,-TRACK_SEP*0.5,"3'→5'",ha='right',va='center',fontsize=8.5*fs,color=C_NEG,fontweight='bold')
    sb_len = 200_000 if span<800_000 else 500_000
    sb_x   = r1-sb_len*1.3
    y_sb   = -(n_neg+1.6)*TRACK_SEP
    ax.hlines(y_sb,sb_x,sb_x+sb_len,colors=C_AXIS,lw=2)
    ax.vlines([sb_x,sb_x+sb_len],y_sb-0.06,y_sb+0.06,colors=C_AXIS,lw=1.5)
    ax.text(sb_x+sb_len/2,y_sb-0.13,f"{sb_len//1000} kb",ha='center',va='top',fontsize=7.5*fs)
    for t in np.linspace(r0,r1,5):
        ax.vlines(t,-0.09,0.09,colors='#999999',lw=0.8)
    ax.text((r0+r1)/2,y_sb-0.65,"chr6 position (Mb)",ha='center',va='top',fontsize=7.5*fs,color='#444444')
    ax.text(r0-span*0.013,(n_pos+1.8)*TRACK_SEP,panel_label,fontsize=13*fs,fontweight='bold',va='bottom')
    ax.set_xlim(r0-span*0.016,r1+span*0.005)
    ax.set_ylim(y_sb-0.85,(n_pos+2.8)*TRACK_SEP)
    ax.axis('off')

# ── Vertical locus panel (new: chromosome backbone runs top→bottom) ─────────────
def draw_panel_v(ax, genes, r0, r1, panel_label, label_fontsize, fs,
                 gap_frac=0.008, min_disp_frac=0.004):
    span   = r1 - r0
    tracks = assign_tracks(genes, gap_frac=gap_frac, region_span=span)
    n_pos  = max((tracks[l] for (_,_,s,l,_) in genes if s=='+'), default=-1)+1
    n_neg  = max((tracks[l] for (_,_,s,l,_) in genes if s=='-'), default=-1)+1
    GW     = GENE_H
    gene_info = {l: (s, tracks[l]) for (_,_,s,l,_) in genes}

    for start, end, strand, label, is_pseudo in genes:
        t     = tracks[label]
        color = gene_color(strand, is_pseudo)
        fc    = color if not is_pseudo else 'white'
        xc    = (t+1)*TRACK_SEP*(1 if strand=='+' else -1)
        dh    = max(end-start, span*min_disp_frac)
        cy    = (start+end)/2
        ys    = cy - dh/2
        ah    = min(span*0.015, dh*0.40)
        bh    = max(dh-ah, dh*0.02)
        xl    = xc - GW/2
        if strand == '+':
            ax.add_patch(plt.Rectangle((xl,ys),GW,bh,fc=fc,ec=color,lw=0.9,zorder=3))
            ax.fill([xl,xl+GW,xc],[ys+bh,ys+bh,ys+dh],fc=fc,ec=color,lw=0.9,zorder=3)
        else:
            ax.add_patch(plt.Rectangle((xl,ys+ah),GW,bh,fc=fc,ec=color,lw=0.9,zorder=3))
            ax.fill([xl,xl+GW,xc],[ys+ah,ys+ah,ys],fc=fc,ec=color,lw=0.9,zorder=3)
        sg = 0.05
        if strand == '+':
            x0 = xc+GW/2; x1 = x0+sg+t*0.10; ha='left'; lx = x1+0.03
        else:
            x0 = xc-GW/2; x1 = x0-sg-t*0.10; ha='right'; lx = x1-0.03
        ax.plot([x0,x1],[cy,cy],color='#bbbbbb',lw=0.5,zorder=2)
        ax.text(lx, cy, label, ha=ha, va='center', fontsize=label_fontsize,
                fontstyle='italic' if is_pseudo else 'normal',
                color='#555555' if is_pseudo else '#111111', clip_on=True)

    ax.vlines(0,r0,r1,colors=C_AXIS,lw=1.8,zorder=2)

    for (yov,l1,l2) in find_overlaps(genes):
        s1,t1 = gene_info[l1]; s2,t2 = gene_info[l2]
        x1c = (t1+1)*TRACK_SEP*(1 if s1=='+' else -1)
        x2c = (t2+1)*TRACK_SEP*(1 if s2=='+' else -1)
        ax.hlines(yov,min(x1c,x2c)-GW/2,max(x1c,x2c)+GW/2,
                  colors=C_OVERLAP,lw=1.5,linestyles='dotted',zorder=6)

    y_top = r0 - span*0.04
    ax.text( TRACK_SEP*0.5, y_top,"5'→3'",ha='center',va='bottom',fontsize=8.5*fs,color=C_POS,fontweight='bold')
    ax.text(-TRACK_SEP*0.5, y_top,"3'→5'",ha='center',va='bottom',fontsize=8.5*fs,color=C_NEG,fontweight='bold')

    sb_len = 200_000 if span<800_000 else 500_000
    x_sb   = (n_pos+1.8)*TRACK_SEP
    y_sb0  = r1 - span*0.12
    ax.vlines(x_sb,y_sb0,y_sb0+sb_len,colors=C_AXIS,lw=2)
    ax.hlines([y_sb0,y_sb0+sb_len],x_sb-0.05,x_sb+0.05,colors=C_AXIS,lw=1.5)
    ax.text(x_sb+0.12,y_sb0+sb_len/2,f"{sb_len//1000} kb",ha='left',va='center',fontsize=7.5*fs)

    for tp in np.linspace(r0,r1,5):
        ax.hlines(tp,-0.08,0.08,colors='#999999',lw=0.8)
    ax.text(-(n_neg+2.5)*TRACK_SEP,(r0+r1)/2,"chr6 position (Mb)",
            ha='center',va='center',fontsize=7.5*fs,color='#444444',rotation=90)
    ax.text(-(n_neg+2.0)*TRACK_SEP,r0-span*0.07,
            panel_label,fontsize=13*fs,fontweight='bold',va='bottom')

    ax.set_xlim(-(n_neg+3.5)*TRACK_SEP,(n_pos+3.5)*TRACK_SEP)
    ax.set_ylim(r1+span*0.05, r0-span*0.12)   # y inverted: top=r0, bottom=r1
    ax.axis('off')

# ── Panel C ────────────────────────────────────────────────────────────────────
GENES_C = ['HLA-A','HLA-DRB1','HLA-B']
PSEUDOGENES_C = {
    'HLA-A':   'pseudogenes: H, J, L, P, V',
    'HLA-DRB1':'pseudogenes: DRB2, DRB6, DRB7, DRB8, DRB9',
    'HLA-B':   'pseudogenes: K, N, T, U, W',
}
BOUNDARIES = {
    'HLA-A':   {'cds_start':22, 'utr3_start':1117},
    'HLA-DRB1':{'cds_start':106,'utr3_start':904},
    'HLA-B':   {'cds_start':21, 'utr3_start':1107},
}
SMOOTH=80; COL_SHARED='#b03a2e'; COL_LINE='#7b241c'; COL_UNIQUE='#d5d8dc'
COL_50='#566573'; COL_UTR='#7fb3d3'; COL_CDS='#2e4053'; COL_BOUND='#1a252f'

def smooth_s(arr,win):
    return pd.Series(arr).rolling(win,center=True,min_periods=5).mean().values

def draw_annotation_c(ax,b,xmax,fs):
    y0,y1=1.04,1.13; bh=y1-y0
    cs,u3=b['cds_start'],b['utr3_start']
    ax.add_patch(mpatches.Rectangle((0,y0),cs,bh,color=COL_UTR,lw=0,zorder=6))
    ax.add_patch(mpatches.Rectangle((cs,y0),u3-cs,bh,color=COL_CDS,lw=0,zorder=6))
    ax.add_patch(mpatches.Rectangle((u3,y0),xmax-u3,bh,color=COL_UTR,lw=0,zorder=6))
    for xb in (cs,u3):
        ax.axvline(xb,color=COL_BOUND,lw=0.9,ls='--',alpha=0.55,zorder=5)
    bm=(y0+y1)/2
    fk=dict(fontsize=7.5*fs,va='center',ha='center',fontweight='bold',color='white',zorder=7)
    if cs>30: ax.text(cs/2,bm,"5'UTR",**fk)
    ax.text((cs+u3)/2,bm,'CDS',**fk)
    if xmax-u3>60: ax.text((u3+xmax)/2,bm,"3'UTR",**fk)
    for xb,lbl in ((cs,f'{cs}'),(u3,f'{u3}')):
        ax.text(xb,-0.07,lbl,ha='center',va='top',fontsize=6.5*fs,color=COL_BOUND,zorder=5)

def draw_panel_c(axes_list, df, fs):
    for i,(gene,ax) in enumerate(zip(GENES_C,axes_list)):
        sub  = df[df['functional_gene']==gene].sort_values('position')
        pos  = sub['position'].values
        flags= sub['shared_90pct'].values.astype(float)
        sm   = smooth_s(flags,SMOOTH)
        b    = BOUNDARIES[gene]; xmax=pos[-1]
        ax.fill_between(pos,0,1, color=COL_UNIQUE,alpha=1.0,lw=0,zorder=1)
        ax.fill_between(pos,0,sm,color=COL_SHARED,alpha=0.85,lw=0,zorder=2)
        ax.plot(pos,sm,color=COL_LINE,lw=1.0,zorder=3)
        ax.axhline(0.5,color=COL_50,lw=0.8,ls='--',zorder=4,alpha=0.7)
        draw_annotation_c(ax,b,xmax,fs)
        pct=100*flags.mean()
        ax.text(0.985,0.70,f'{pct:.0f}% non-unique',transform=ax.transAxes,
                ha='right',va='top',fontsize=9*fs,color=COL_LINE,fontweight='bold')
        ax.set_title(f'{gene}   ({PSEUDOGENES_C[gene]})',loc='left',
                     fontsize=8*fs,fontweight='bold',pad=14)
        ax.set_xlim(0, xmax + xmax * 0.012); ax.set_ylim(-0.04,1.16)
        ax.set_yticks([0,0.25,0.5,0.75,1.0])
        ax.set_yticklabels(['0%','25%','50%','75%','100%'],fontsize=8*fs)
        ax.set_ylabel('shared ≥90%\nwith pseudogene',fontsize=8*fs,labelpad=4)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#aab7b8'); ax.spines['bottom'].set_color('#aab7b8')
        ax.tick_params(axis='both',colors='#566573',labelsize=8*fs)
        if i<len(GENES_C)-1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('Position in MANE Select transcript (bp, 100 bp window start)',
                          fontsize=9*fs, labelpad=12*fs)
    axes_list[0].text(-0.12,1.22,'C',transform=axes_list[0].transAxes,
                      fontsize=13*fs,fontweight='bold',va='top')

def c_legend(ax, fs):
    handles=[
        mpatches.Patch(color=COL_SHARED,alpha=0.85,label='≥90% identical to ≥1 pseudogene (100 bp window)'),
        mpatches.Patch(color=COL_UNIQUE,label='unique region'),
        mpatches.Patch(color=COL_CDS,label='CDS'),
        mpatches.Patch(color=COL_UTR,label="5'/3' UTR"),
        plt.Line2D([0],[0],color=COL_50,lw=0.9,ls='--',label='50%'),
    ]
    ax.legend(handles=handles,loc='upper center',ncol=3,fontsize=8*fs,
              frameon=False,bbox_to_anchor=(0.5,-0.38))

df_shared = pd.read_csv(DIR_MAFFT/'shared_positions.tsv', sep='\t')

# ══════════════════════════════════════════════════════════════════════════════
# Figure_ABC: A+B stacked left | C right
# ══════════════════════════════════════════════════════════════════════════════
fig1 = plt.figure(figsize=(36, 20))
gs1  = gridspec.GridSpec(1, 2, figure=fig1,
                         width_ratios=[1.15, 1.0],
                         wspace=0.28)          # gap between A/B column and C

left_gs1 = gridspec.GridSpecFromSubplotSpec(
    2, 1, subplot_spec=gs1[0], hspace=0.30, height_ratios=[1.4, 1.0])
ax_a1 = fig1.add_subplot(left_gs1[0])
ax_b1 = fig1.add_subplot(left_gs1[1])

draw_panel_h(ax_a1, PANEL_A_GENES, 32_750_000, 33_400_000, 'A',
             6.5*FS_AB, FS_AB, gap_frac=0.003, min_disp_frac=0.004)
draw_panel_h(ax_b1, PANEL_B_GENES, 29_650_000, 31_600_000, 'B',
             7.5*FS_AB, FS_AB, gap_frac=0.008, min_disp_frac=0.003)

# B legend: 2 rows, 1.25x font, lowered to align with C legend
ax_b1.legend(handles=locus_legend_handles(), loc='lower center', ncol=2,
             fontsize=6.25*FS_AB, framealpha=0.9,
             handlelength=0.9, handleheight=0.6,
             borderpad=0.4, labelspacing=0.3,
             bbox_to_anchor=(0.5, -0.26))

right_gs1 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs1[1], hspace=0.48)
axes_c1   = [fig1.add_subplot(right_gs1[i]) for i in range(3)]

draw_panel_c(axes_c1, df_shared, fs=FS_C)
c_legend(axes_c1[-1], FS_C)

out1 = Path(__file__).parent / 'Figure_ABC'
for fmt in ('pdf','png'):
    fig1.savefig(f'{out1}.{fmt}', dpi=200, bbox_inches='tight')
    print(f'Saved → {out1}.{fmt}')
plt.close(fig1)

# ══════════════════════════════════════════════════════════════════════════════
# Figure_ABC_vchrom: vertical chr lines  [A | B | C]
# ══════════════════════════════════════════════════════════════════════════════
fig2 = plt.figure(figsize=(34, 20))
gs2  = gridspec.GridSpec(1, 3, figure=fig2,
                         width_ratios=[0.6, 0.55, 1.0],
                         wspace=0.18)

ax_a2 = fig2.add_subplot(gs2[0])
ax_b2 = fig2.add_subplot(gs2[1])
cgs2  = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs2[2], hspace=0.55)
axes_c2 = [fig2.add_subplot(cgs2[i]) for i in range(3)]

draw_panel_v(ax_a2, PANEL_A_GENES, 32_750_000, 33_400_000, 'A',
             6.5*FS_AB, FS_AB, gap_frac=0.003, min_disp_frac=0.004)
draw_panel_v(ax_b2, PANEL_B_GENES, 29_650_000, 31_600_000, 'B',
             7.5*FS_AB, FS_AB, gap_frac=0.008, min_disp_frac=0.003)

fig2.legend(handles=locus_legend_handles(), loc='lower left',
            ncol=4, fontsize=8*FS_AB, framealpha=0.9,
            bbox_to_anchor=(0.01, 0.01))

draw_panel_c(axes_c2, df_shared, fs=FS)
c_legend(axes_c2[-1], FS)

out2 = Path(__file__).parent / 'Figure_ABC_vchrom'
for fmt in ('pdf','png'):
    fig2.savefig(f'{out2}.{fmt}', dpi=200, bbox_inches='tight')
    print(f'Saved → {out2}.{fmt}')
plt.close(fig2)
