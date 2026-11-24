#!/usr/bin/env python3
"""
Overview figure of mappability phenomena, from existing results only (no new runs).
  A  funique per reference, all MANE transcripts (ECDF), 100-bp PE simulation
  B  APM genes: funique per reference (heat map)
  C  kallisto TPM on the dedup-CDS index / truth (uniform reads), all transcripts
  D  number of CDS isoforms vs kallisto ratio
  E  HLA genes: % windows shared with pseudogenes vs read length (HLApseudo_MAFFT)
  F  RNU snRNA families: funique (SE 100 bp)
Run from mappability-correction/:  python3 Collesano2026/scripts/plot_phenomena.py
Outputs: Collesano2026/results/figures/phenomena_overview.{png,pdf}, results/tables/phenomena_*.tsv
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

R = Path("mappability/results")
OUT_F, OUT_T = Path("Collesano2026/results/figures"), Path("Collesano2026/results/tables")
OUT_F.mkdir(parents=True, exist_ok=True); OUT_T.mkdir(parents=True, exist_ok=True)

sym = {}
for ln in open("Ref/MANE.GRCh38.v1.5.ensembl_rna.fna"):
    if ln[0] == ">":
        sym[ln[1:].split()[0]] = re.search(r"gene_symbol:(\S+)", ln).group(1)

REFS = {"genome + alt (RefSeq p14)": R / "transcript_uniqueness_factors_genomic_L100bp.tsv",
        "dedup CDS": R / "transcript_uniqueness_factors_dedup_cds_RNA_PE.tsv",
        "MANE only": R / "transcript_uniqueness_factors_MANE_RNA_PE.tsv",
        "Ensembl primary (used)": Path("Collesano2026/resources/funique/transcript_uniqueness_factors_L100bp.tsv")}
COL = {"genome + alt (RefSeq p14)": "#d62728", "dedup CDS": "#ff7f0e", "MANE only": "#7f7f7f",
       "Ensembl primary (used)": "#1f77b4"}
uf = {}
for k, f in REFS.items():
    t = pd.read_csv(f, sep="\t")
    uf[k] = t.groupby(t.transcript_id.map(sym)).uniqueness_factor.first()
uf = pd.DataFrame(uf)

APM = ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1",
       "HLA-DQB1", "HLA-DRB1", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "MICA", "MICB", "B2M"]

k = pd.read_csv("mappability/kallisto_results/kallisto_inflation.tsv", sep="\t")
k = k[(k.read_len_bp == 100) & (k["mode"] == "PE") & (k.read_source == "MANE_RNA")].copy()
k = k[k.tpm_mane > 0]                                   # truth index must see the transcript
k["ratio"] = k.tpm_cds / k.tpm_mane
iso = pd.read_csv("mappability/smallPhenomena/01_isoform_redundancy/isoform_table.tsv", sep="\t")
k = k.merge(iso[["gene_name", "n_isoforms"]], left_on="gene_symbol", right_on="gene_name", how="left")
sh = pd.read_csv("HLApseudo_MAFFT/shared_fraction.tsv", sep="\t")
rnu = pd.read_csv("RNU_snRNA/RNU_UF_summary.tsv", sep="\t")

# tables
summ = pd.DataFrame({r: {"n": uf[r].notna().sum(), "f=0": (uf[r] == 0).sum(), "f<0.9": (uf[r] < .9).sum(),
                         "median": uf[r].median()} for r in uf}).T
summ.to_csv(OUT_T / "phenomena_funique_by_reference.tsv", sep="\t")
uf.reindex(APM).to_csv(OUT_T / "phenomena_APM_funique_by_reference.tsv", sep="\t", float_format="%.3f")
kz = k.assign(cls=pd.cut(k.ratio, [-1, 1e-6, .9, 1.1, np.inf], labels=["zero", "lower >10%", "within 10%", "higher >10%"]))
kz.cls.value_counts().to_csv(OUT_T / "phenomena_kallisto_dedupCDS_classes.tsv", sep="\t")

fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))
a = ax[0, 0]
for r in uf:
    x = np.sort(uf[r].dropna().values)
    a.step(x, np.arange(1, len(x) + 1) / len(x), color=COL[r], label=f"{r} (f=0: {int(summ.loc[r, 'f=0'])})")
a.set_xlabel("funique (MAPQ-255 fraction)"); a.set_ylabel("cumulative fraction of MANE transcripts")
a.set_title("A  same reads, four references", loc="left", fontsize=10); a.legend(fontsize=7, frameon=False)

a = ax[0, 1]
m = uf.reindex(APM)
im = a.imshow(m.values, cmap="RdYlBu", vmin=0, vmax=1, aspect="auto")
a.set_yticks(range(len(APM))); a.set_yticklabels(APM, fontsize=7)
a.set_xticks(range(m.shape[1])); a.set_xticklabels([c.replace(" (", "\n(") for c in m.columns], fontsize=7)
for i in range(m.shape[0]):
    for j in range(m.shape[1]):
        v = m.values[i, j]
        if pd.notna(v):
            a.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6)
a.set_title("B  APM genes: alt contigs give zero", loc="left", fontsize=10); plt.colorbar(im, ax=a, fraction=.04)

a = ax[0, 2]
x = np.log2(k.ratio.clip(lower=2 ** -8))
a.hist(x, bins=80, color="#999999")
for g in ["TAPBP", "TAP2", "HLA-A", "B2M", "HLA-G"]:
    v = k.loc[k.gene_symbol == g, "ratio"]
    if len(v):
        a.axvline(np.log2(max(v.iloc[0], 2 ** -8)), color="#d62728", lw=.8)
        a.text(np.log2(max(v.iloc[0], 2 ** -8)), a.get_ylim()[1] * .9, g, rotation=90, fontsize=7, va="top")
a.set_yscale("log"); a.set_xlabel("log2(kallisto TPM, dedup-CDS index / truth)")
a.set_title(f"C  EM gives fractions (zero {100 * (k.ratio == 0).mean():.0f} %)", loc="left", fontsize=10)

a = ax[1, 0]
d = k.dropna(subset=["n_isoforms"]).copy()
d["iso"] = d.n_isoforms.clip(upper=10).astype(int)
grp = [d.loc[d.iso == i, "ratio"].values for i in range(1, 11)]
a.boxplot(grp, positions=range(1, 11), showfliers=False)
a.set_xticks(range(1, 11)); a.set_xticklabels([str(i) for i in range(1, 10)] + ["≥10"])
a.set_xlabel("CDS isoforms per gene"); a.set_ylabel("kallisto ratio (dedup CDS / truth)")
a.set_title("D  isoforms split the signal", loc="left", fontsize=10)

a = ax[1, 1]
s = sh[sh.threshold_pct == 90]
for g in ["HLA-A", "HLA-G", "HLA-B", "HLA-C", "HLA-E", "HLA-DRB1"]:
    t = s[s.functional_gene == g]
    a.plot(t.window_bp, 100 * t.frac_shared, "o-", label=g)
a.set_xlabel("read length (window, bp)"); a.set_ylabel("% windows ≥ 90 % identical to a pseudogene")
a.set_title("E  HLA vs pseudogenes (lower bound)", loc="left", fontsize=10); a.legend(fontsize=7, frameon=False)

a = ax[1, 2]
r = rnu.copy(); r["SE_100bp"] = pd.to_numeric(r.SE_100bp, errors="coerce")
a.barh(r.gene_symbol, r.SE_100bp, color=["#d62728" if v == 0 else "#1f77b4" for v in r.SE_100bp])
a.invert_yaxis(); a.tick_params(axis="y", labelsize=6); a.set_xlabel("funique (MANE, SE 100 bp)")
a.set_title("F  RNU snRNA copies: identical loci = 0", loc="left", fontsize=10)

plt.tight_layout()
for ext in ("png", "pdf"):
    plt.savefig(OUT_F / f"phenomena_overview.{ext}", dpi=150)
print(summ.to_string()); print(kz.cls.value_counts().to_string())
