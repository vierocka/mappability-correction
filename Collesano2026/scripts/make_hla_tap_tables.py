#!/usr/bin/env python3
"""
Regenerate the HLA and TAP1/TAP2 summary tables from corrected_counts/.

Run AFTER apply_uniqueness_and_average.py.  Reads only its outputs, so the
numbers here can never drift from the corrected count files.

Outputs (schemas unchanged from the existing files, with `f_unique` appended
to the HLA table because that column is what the 2026-10-06 recompute changed):

  HLA_corrected_mean_values.tsv        per condition x HLA allele
  HLA_per_sample_trace.tsv             per sample x HLA allele
  HLA_gene_level_summary.tsv           per condition x gene, alleles summed
  TAPs/TAP1_TAP2_condition_mean_trace.tsv
  TAPs/TAP1_TAP2_per_sample_trace.tsv

Allele rows are summed in HLA_gene_level_summary.tsv because a heterozygous
sample's two alleles are two rows of the SAME gene: comparing a single allele
row across cell lines compares different quantities (JY is homozygous and has
one row carrying the whole gene's signal, HCC1937/HCT116 split theirs in two).
"""

from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent

# --- PIPELINE HOOK (2026-11-22): base dir + config from Snakemake ---------
import json as _json, os as _os
BASE = Path(_os.environ.get("PIPE_BASE", BASE))
_CFG = _json.load(open(_os.environ["PIPE_CONFIG"])) if _os.environ.get("PIPE_CONFIG") else {}
# -------------------------------------------------------------------------

HLA_BASE = {"ENST00000376809": "HLA-A", "ENST00000412585": "HLA-B",
            "ENST00000376228": "HLA-C"}
TAPS = {"ENST00000354258.5": "TAP1", "ENST00000374897.4": "TAP2"}

CONDITIONS = [                      # cell, condition, [samples]
    ("HCC1937", "74_75",    ["SRR31477874", "SRR31477875"]),   # SRR925709: older batch, not mixed
    ("HCT116",  "41_42",    ["SRR3644141", "SRR3644142"]),
    ("JY",      "364065",   ["SRR364065"]),
]


# PIPELINE HOOK
CONDITIONS = [tuple(c) for c in _CFG.get("CONDITIONS", CONDITIONS)]


def gene_of(tid):
    base = tid.split("_")[0].split(".")[0]
    return HLA_BASE.get(base)


def mean_path(cell, cond):
    return BASE / cell / "corrected_counts" / f"{cell}_corrected_mean_{cond}_mane_fc_corrected.tsv"


def sample_path(cell, srr):
    return BASE / cell / "corrected_counts" / f"{srr}_mane_fc_corrected.tsv"


# f_unique is a property of the cell line, identical across its samples
fu = {}
for cell, cond, samples in CONDITIONS:
    d = pd.read_csv(sample_path(cell, samples[0]), sep="\t")
    fu[cell] = d.set_index("gene_id")["f_unique"].to_dict()

hla_rows, hla_samp, gene_rows, tap_cond, tap_samp = [], [], [], [], []

for cell, cond, samples in CONDITIONS:
    m = pd.read_csv(mean_path(cell, cond), sep="\t")
    rel = mean_path(cell, cond).relative_to(BASE).as_posix()

    # ── HLA, per condition x allele ─────────────────────────────────────────
    h = m[m.gene_id.map(gene_of).notna()].copy()
    h["gene"] = h.gene_id.map(gene_of)
    for _, r in h.iterrows():
        hla_rows.append({
            "PATH:original_MANE_protein_ID_PredictedHLAallele": f"{rel}:{r.gene_id}",
            "original_length": r.length,
            "mean_counts": r.mean_counts,
            "mean_counts_corrected": r.mean_counts_corrected,
            "mean_TPM": r.mean_TPM,
            "mean_TPM_corrected": r.mean_TPM_corrected,
            "f_unique": fu[cell].get(r.gene_id),
            "_gene": r.gene, "_cell": cell, "_cond": cond,
        })

    # ── HLA, per condition x gene (alleles summed) ──────────────────────────
    for gene, g in h.groupby("gene"):
        gene_rows.append({
            "cell_line": cell, "condition": cond, "gene": gene,
            "n_alleles": len(g),
            "alleles": ",".join(x.split("_")[1] if "_" in x else "homozygous"
                                for x in g.gene_id),
            "mean_counts": g.mean_counts.sum(),
            "mean_counts_corrected": g.mean_counts_corrected.sum(),
            "mean_TPM": g.mean_TPM.sum(),
            "mean_TPM_corrected": g.mean_TPM_corrected.sum(),
        })

    # ── TAP1/TAP2, per condition ────────────────────────────────────────────
    for tid, gene in TAPS.items():
        r = m[m.gene_id == tid]
        if len(r):
            r = r.iloc[0]
            tap_cond.append({"cell_line": cell, "condition": cond, "gene": gene,
                             "transcript_id": tid,
                             "mean_counts": r.mean_counts,
                             "mean_counts_corrected": r.mean_counts_corrected,
                             "mean_TPM": r.mean_TPM,
                             "mean_TPM_corrected": r.mean_TPM_corrected})

    # ── per-sample traces ───────────────────────────────────────────────────
    for srr in samples:
        s = pd.read_csv(sample_path(cell, srr), sep="\t")
        for _, r in s[s.gene_id.map(gene_of).notna()].iterrows():
            hla_samp.append({"cell_line": cell, "sample": srr,
                             "gene": gene_of(r.gene_id), "transcript_id": r.gene_id,
                             "length": r.length, "raw_counts": r.counts,
                             "f_unique_applied": r.f_unique,
                             "counts_corrected": r.counts_corrected,
                             "TPM_raw": r.TPM, "TPM_corrected": r.TPM_corrected})
        for tid, gene in TAPS.items():
            r = s[s.gene_id == tid]
            if len(r):
                r = r.iloc[0]
                tap_samp.append({"cell_line": cell, "sample": srr, "gene": gene,
                                 "transcript_id": tid, "length": r.length,
                                 "raw_counts": r.counts,
                                 "f_unique_applied": r.f_unique,
                                 "counts_corrected": r.counts_corrected,
                                 "TPM_raw": r.TPM, "TPM_corrected": r.TPM_corrected})

# ── write, preserving the original row order (gene, then cell/condition) ────
hla = pd.DataFrame(hla_rows).sort_values(["_gene", "_cell", "_cond",
                                          "PATH:original_MANE_protein_ID_PredictedHLAallele"])
hla.drop(columns=["_gene", "_cell", "_cond"]).to_csv(
    BASE / "HLA_corrected_mean_values.tsv", sep="\t", index=False, float_format="%.6f")

pd.DataFrame(hla_samp).to_csv(BASE / "HLA_per_sample_trace.tsv", sep="\t",
                              index=False, float_format="%.4f")
pd.DataFrame(gene_rows).sort_values(["gene", "cell_line", "condition"]).to_csv(
    BASE / "HLA_gene_level_summary.tsv", sep="\t", index=False, float_format="%.4f")
(BASE / "TAPs").mkdir(exist_ok=True)
pd.DataFrame(tap_cond).to_csv(BASE / "TAPs/TAP1_TAP2_condition_mean_trace.tsv",
                              sep="\t", index=False, float_format="%.4f")
pd.DataFrame(tap_samp).to_csv(BASE / "TAPs/TAP1_TAP2_per_sample_trace.tsv",
                              sep="\t", index=False, float_format="%.4f")

print("Wrote HLA_corrected_mean_values.tsv, HLA_per_sample_trace.tsv, "
      "HLA_gene_level_summary.tsv, TAPs/TAP1_TAP2_{condition_mean,per_sample}_trace.tsv\n")
print(pd.DataFrame(gene_rows).sort_values(["gene", "cell_line", "condition"])
      .to_string(index=False))
