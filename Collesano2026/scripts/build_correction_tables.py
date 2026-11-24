#!/usr/bin/env python3
"""
Build ONE uniqueness-correction table per read length, covering every MANE
Select transcript.

Scope: basic MANE only -- one transcript per gene, HLA-A/B/C included (the
GRCh38 reference allele). Identical for every sample; nothing cell-line- or
allele-specific. Heterozygosity is handled downstream, not here.

Source: the primary-assembly base run (HPC 2026-08-12), 19,373 transcripts,
plus explicit rows for the 64 MANE transcripts whose locus is not on a primary
chromosome and for which no measurement exists. 19,373 + 64 = 19,437 = all of
MANE Select.

Output: correction_tables/uniqueness_correction_L{75,100}bp.tsv
"""

import re
from pathlib import Path

import pandas as pd

MANE = Path(__file__).resolve().parents[3] / "MANE"
AUG = MANE / "primary_assembly_Aug12"
OUT = MANE / "correction_tables"

# --- PIPELINE HOOK (2026-11-22): any read length, paths from Snakemake ------
import os as _os
MANE = Path(_os.environ.get("PIPE_MANE_DIR", MANE))          # MANE rna.fna + genomic.gff
AUG = Path(_os.environ.get("PIPE_FACTORS_DIR", AUG))         # transcript_uniqueness_factors_L{L}bp.tsv
OUT = Path(_os.environ.get("PIPE_OUT_DIR", OUT))
LENGTHS = [int(x) for x in _os.environ.get("PIPE_READ_LENS", "75,100").split(",")]
OUT.mkdir(parents=True, exist_ok=True)

# ── MANE annotation ─────────────────────────────────────────────────────────
symbol = {}
for line in open(MANE / "MANE.GRCh38.v1.5.ensembl_rna.fna"):
    if line.startswith(">"):
        tid = line[1:].split()[0]
        symbol[tid] = re.search(r"gene_symbol:(\S+)", line).group(1)

primary = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]}
non_primary = {}
for line in open(MANE / "MANE.GRCh38.v1.5.ensembl_genomic.gff"):
    f = line.split("\t")
    if len(f) > 8 and f[2] == "transcript" and f[0] not in primary:
        tid = re.search(r"transcript_id=([^;\n]+)", f[8]).group(1)
        non_primary[tid] = ("alt_haplotype" if f[0].endswith("_alt") else "fix_patch")

COLS = ["transcript_id", "gene_symbol", "read_len", "n_positions",
        "n_unique_back", "n_multi_back", "n_unmapped", "uniqueness_factor", "status"]

for L in LENGTHS:
    base = pd.read_csv(AUG / f"transcript_uniqueness_factors_L{L}bp.tsv", sep="\t")
    base["gene_symbol"] = base.transcript_id.map(symbol)
    base["read_len"] = L
    base["status"] = base.uniqueness_factor.map(lambda v: "zero" if v == 0 else "ok")

    missing = [t for t in symbol if t not in set(base.transcript_id)]
    exc = pd.DataFrame([{
        "transcript_id": t, "gene_symbol": symbol[t], "read_len": L,
        "n_positions": pd.NA, "n_unique_back": pd.NA, "n_multi_back": pd.NA,
        "n_unmapped": pd.NA, "uniqueness_factor": pd.NA,
        "status": f"excluded_{non_primary[t]}",
    } for t in missing])

    tab = pd.concat([base[COLS], exc[COLS]], ignore_index=True)
    tab = tab.sort_values(["gene_symbol", "transcript_id"], kind="stable")

    assert len(tab) == 19437, f"expected 19,437 MANE transcripts, got {len(tab)}"
    assert not tab.transcript_id.duplicated().any(), "transcript_id not unique"

    path = OUT / f"uniqueness_correction_L{L}bp.tsv"
    tab.to_csv(path, sep="\t", index=False)
    n = tab.status.value_counts()
    print(f"{path.name}: {len(tab):,} rows  "
          f"(ok {n.get('ok', 0):,}, zero {n.get('zero', 0)}, "
          f"no factor {n.get('excluded_fix_patch', 0) + n.get('excluded_alt_haplotype', 0)})")
    hla = tab[tab.gene_symbol.isin(["HLA-A", "HLA-B", "HLA-C"])]
    print("   " + "  ".join(f"{r.gene_symbol} {r.uniqueness_factor:.4f}"
                            for _, r in hla.iterrows()))
