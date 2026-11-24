#!/usr/bin/env python3
"""
Rebuild the peptide-level EXPRESSION tables from the refreshed corrected_counts/.
Run AFTER apply_uniqueness_and_average.py.

Scope: expression only.  P_chop and the GEM energies (en_0201/en_0702) belong
to the chopping/energy side of the project and are deliberately NOT produced
or carried here -- including anything derived from them (the old TPM_prod
column).  Both remain available to join back:

  P_chop   -> {CELL}_9mer_pchop.csv            on (transcript_id, peptide)
  en_*     -> peptide_abundance_{CELL}_PRE<date>.csv.bak   on (peptide)

NOTE for whoever owns TPM_prod: the previous files were inconsistent between
cell lines.  JY's TPM_prod was P_chop x RAW mean_TPM; HCC1937's and HCT116's
were P_chop x mean_TPM_corrected (verified against the backups: ratio 1.000000
for JY, 0.9994/0.9996 for the other two).  Cross-cell-line comparisons on that
column were mixing raw and corrected values, differing by up to ~1000x on
transcripts floored at f_unique = 0.001.  Recompute it from one definition.

Outputs:
  master_file_abundance_{CELL}.csv  transcript_id, gene_symbol, peptide,
                                    mean_TPM, mean_TPM_corrected
  peptide_abundance_{CELL}.csv      peptide, mean_TPM, mean_TPM_corrected
                                    (summed over transcripts carrying it)

Peptide-level rows keep the original filters: mean_TPM_corrected > 0, and
peptides containing U or X dropped.
"""

import shutil
from datetime import date
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent

# --- PIPELINE HOOK (2026-11-22): base dir + config from Snakemake ---------
import json as _json, os as _os
BASE = Path(_os.environ.get("PIPE_BASE", BASE))
_CFG = _json.load(open(_os.environ["PIPE_CONFIG"])) if _os.environ.get("PIPE_CONFIG") else {}
# -------------------------------------------------------------------------
CHUNK = 2_000_000
STAMP = date.today().isoformat()

EXPR = {
    "JY":      "JY_corrected_mean_364065_mane_fc_corrected.tsv",
    "HCT116":  "HCT116_corrected_mean_41_42_mane_fc_corrected.tsv",
    "HCC1937": "HCC1937_corrected_mean_74_75_mane_fc_corrected.tsv",
}

# PIPELINE HOOK: expr files + 9-mer source per cell line
EXPR = _CFG.get("EXPR", EXPR)
PCHOP = _CFG.get("PCHOP", {})   # cell -> transcript x 9-mer CSV (one row per position)
MIN_TPM = float(_CFG.get("min_transcript_tpm", 0))   # transcript TPM_corrected cut-off (pipeline default 1)

for cell, expr_name in EXPR.items():
    print(f"\n{'=' * 60}\n{cell}\n{'=' * 60}")
    cdir = BASE / cell
    expr = pd.read_csv(cdir / "corrected_counts" / expr_name, sep="\t")
    expr = expr.rename(columns={"gene_id": "transcript_id"})
    tpm = expr.set_index("transcript_id")["mean_TPM"].to_dict()
    tpm_c = expr.set_index("transcript_id")["mean_TPM_corrected"].to_dict()

    master = cdir / f"master_file_abundance_{cell}.csv"
    for path in (master, cdir / f"peptide_abundance_{cell}.csv"):
        bak = path.with_name(path.stem + f"_PRE{STAMP}.csv.bak")
        if path.exists() and not bak.exists():
            shutil.copy2(path, bak)
            print(f"  backed up {path.name} -> {bak.name}")

    tmp = master.with_suffix(".csv.tmp")
    n_in = n_out = 0
    parts = []
    with open(tmp, "w") as out:
        for i, ch in enumerate(pd.read_csv(PCHOP.get(cell, cdir / f"{cell}_9mer_pchop.csv"),
                                           usecols=["transcript_id", "gene_symbol", "peptide"],
                                           chunksize=CHUNK)):
            n_in += len(ch)
            ch["mean_TPM"] = ch.transcript_id.map(tpm)
            ch["mean_TPM_corrected"] = ch.transcript_id.map(tpm_c)
            ch = ch[ch.mean_TPM.notna()]
            ch.to_csv(out, index=False, header=(i == 0))
            n_out += len(ch)
            keep = (ch.mean_TPM_corrected >= MIN_TPM) if MIN_TPM > 0 else (ch.mean_TPM_corrected > 0)   # PIPELINE HOOK
            parts.append(ch.loc[keep,
                                ["peptide", "mean_TPM", "mean_TPM_corrected"]])
    tmp.replace(master)
    print(f"  master_file_abundance_{cell}.csv: {n_out:,} / {n_in:,} rows")

    pep = (pd.concat(parts, ignore_index=True)
             .groupby("peptide", as_index=False)[["mean_TPM", "mean_TPM_corrected"]].sum())
    del parts
    pep = pep[~pep.peptide.astype(str).str.contains("U|X", regex=True)].reset_index(drop=True)

    out_pep = cdir / f"peptide_abundance_{cell}.csv"
    bak = out_pep.with_name(out_pep.stem + f"_PRE{STAMP}.csv.bak")
    if bak.exists():
        old = set(pd.read_csv(bak, usecols=["peptide"]).peptide)
        if old != set(pep.peptide):
            raise SystemExit(f"{cell}: peptide set changed -- investigate before overwriting")
        print(f"  peptide set identical to backup ({len(old):,})")
    pep.to_csv(out_pep, index=False)
    print(f"  -> {out_pep.name}  [peptide, mean_TPM, mean_TPM_corrected]  {len(pep):,} rows")

print("\nDone.  P_chop, en_*, TPM_prod and MS_total_expression_* left to their owners.")
