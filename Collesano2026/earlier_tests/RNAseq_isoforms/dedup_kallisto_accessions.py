#!/usr/bin/env python3
import re
import sys
import pandas as pd
from pathlib import Path

ACC_RE = re.compile(r'(NP_\d+\.\d+|XP_\d+\.\d+)')

def dedup_kallisto_by_accession(abundance_tsv: str) -> pd.DataFrame:
    df = pd.read_csv(abundance_tsv, sep="\t")
    # kallisto columns are typically: target_id length eff_length est_counts tpm
    if "target_id" not in df.columns:
        raise ValueError("Expected a kallisto abundance.tsv with a 'target_id' column.")

    def extract_acc(tid: str):
        m = ACC_RE.search(tid)
        return m.group(1) if m else None

    df["accession"] = df["target_id"].map(extract_acc)
    df = df[df["accession"].notna()].copy()

    # Max across patch/ALT duplicates for the same accession
    out = (
        df.groupby("accession", as_index=False)
          .agg(
              est_counts_max=("est_counts", "max"),
              tpm_max=("tpm", "max"),
              n_targets=("target_id", "count")
          )
          .sort_values("est_counts_max", ascending=False)
    )
    return out

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: dedup_kallisto_accessions.py <abundance.tsv> [out.tsv]", file=sys.stderr)
        sys.exit(1)

    in_tsv = sys.argv[1]
    out_tsv = sys.argv[2] if len(sys.argv) >= 3 else None

    dedup = dedup_kallisto_by_accession(in_tsv)

    if out_tsv:
        dedup.to_csv(out_tsv, sep="\t", index=False)
    else:
        print(dedup.to_csv(sep="\t", index=False), end="")

