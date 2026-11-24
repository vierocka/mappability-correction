#!/usr/bin/env python3
"""Does removing non-coding transcripts and rescaling make poly(A) and Ribo-Zero libraries comparable?

The three cell lines were not sequenced with the same chemistry: JY and HCC1937 are poly(A)-selected,
HCT116 is total RNA (Ribo-Zero). TPM is compositional, so a library whose mass sits in non-coding RNA
puts every protein-coding gene on a different scale. The pipeline removes 70 non-coding MANE transcripts
and renormalises. This script asks whether that is enough.

Three diagnostics, all independent of which genes a given cell line happens to express strongly:

  1. non-coding share of TPM, before removal (how much rescaling is being done)
  2. concentration of the protein-coding distribution after rescaling (top-N share, median, IQR)
  3. replication-dependent histones as a chemistry marker. They are the one protein-coding family
     with no poly(A) tail, so they are depleted by poly(A) selection and retained by Ribo-Zero.
     If their share still separates the libraries after rescaling, residual chemistry bias remains
     inside the protein-coding set and cross-line net sums are not on a common scale.

usage: compare_polya_ribozero.py --counts-root DIR [--out DIR]
       DIR holds the published per-run featureCounts tables (see RUNS below).
"""
import argparse
import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
NONCODING = HERE.parent.parent / "resources" / "non_coding_MANE_transcripts_no_protein.tsv"

# run label -> (path relative to --counts-root, chemistry)
RUNS = {
    "JY":           ("JY/SRR364056/MANE/STARfc/SRR364065_JY_mane_fc.counts.csv", "poly(A)"),
    "HCC1937_r1":   ("HCC1937/SRR31477874/MANE/STARfc/SRR31477874_HCC1937_mane_fc.counts.csv", "poly(A)"),
    "HCC1937_r2":   ("HCC1937/SRR31477875/MANE/STARfc/SRR31477875_HCC1937_mane_fc.counts.csv", "poly(A)"),
    "HCT116_r1":    ("HCT116/SRR3644141/MANE/STARfc/SRR3644141_HCT116_mane_fc.counts.csv", "Ribo-Zero"),
    "HCT116_r2":    ("HCT116/SRR3644142/MANE/STARfc/SRR3644142_HCT116_mane_fc.counts.csv", "Ribo-Zero"),
}

# replication-dependent histone symbols: H1-1..H1-10, H2AC*, H2BC*, H3C*, H4C* (HGNC 2019 nomenclature).
# H1-0/H1-10 (HFs), H2AZ/H2AX/H3-3/H3-5 and macroH2A are replication-INDEPENDENT and polyadenylated:
# they are excluded, since they carry a poly(A) tail and are not a chemistry marker.
RD_HISTONE = re.compile(r"^(H1-[1-5]$|H2AC\d|H2BC\d|H3C\d|H4C\d)")
REPL_INDEP = {"H1-0", "H1-10", "H2AZ1", "H2AZ2", "H2AX", "H3-3A", "H3-3B", "H3-4", "H3-5", "MACROH2A1"}


def read_counts(path):
    d = pd.read_csv(path, sep="\t", comment="#")
    d = d.rename(columns={d.columns[0]: "transcript_id", d.columns[5]: "length"})
    d["count"] = pd.to_numeric(d[d.columns[-1]], errors="coerce")
    return d[["transcript_id", "length", "count"]].dropna()


def tpm(d):
    r = d["count"] / (d["length"] / 1000.0)
    return 1e6 * r / r.sum()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--counts-root", required=True)
    p.add_argument("--gene-map", help="TSV with transcript_id and gene_name; "
                                      "default: the non-coding list plus MANE symbols if available")
    p.add_argument("--out", default=str(HERE))
    a = p.parse_args()

    root = Path(a.counts_root)
    nc = pd.read_csv(NONCODING, sep="\t")
    nc_ids = set(nc["transcript_id"])

    gene_name = {}
    if a.gene_map:
        g = pd.read_csv(a.gene_map, sep="\t")
        gene_name = dict(zip(g["transcript_id"], g["gene_name"]))

    rows, hist_rows = [], []
    for label, (rel, chem) in RUNS.items():
        path = root / rel
        if not path.exists():
            print(f"missing, skipped: {path}")
            continue
        d = read_counts(path)
        d["tpm_raw"] = tpm(d)

        is_nc = d.transcript_id.isin(nc_ids)
        nc_share = d.loc[is_nc, "tpm_raw"].sum() / 1e6

        pc = d[~is_nc].copy()
        pc["tpm"] = tpm(pc)                       # rescaled over protein-coding only
        s = pc["tpm"].sort_values(ascending=False)
        expressed = s[s > 0]

        row = {
            "library": label, "chemistry": chem,
            "n_transcripts": len(d),
            "noncoding_pct_of_raw_TPM": round(100 * nc_share, 2),
            "rescale_factor": round(1 / (1 - nc_share), 3),
            "top10_pct_after": round(100 * s.head(10).sum() / 1e6, 2),
            "top100_pct_after": round(100 * s.head(100).sum() / 1e6, 2),
            "median_TPM_expressed": round(expressed.median(), 3),
            "pct_above_TPM1": round(100 * (s >= 1).mean(), 1),
        }

        if gene_name:
            pc["gene"] = pc.transcript_id.map(gene_name)
            sym = pc["gene"].fillna("")
            rd = sym.str.match(RD_HISTONE) & ~sym.isin(REPL_INDEP)
            row["n_RD_histone"] = int(rd.sum())
            row["RD_histone_pct_of_PC_TPM"] = round(100 * pc.loc[rd, "tpm"].sum() / 1e6, 3)
            hist_rows.append(pc.loc[rd, ["gene", "tpm"]].assign(library=label))
        rows.append(row)

    out = pd.DataFrame(rows)
    dest = Path(a.out)
    dest.mkdir(parents=True, exist_ok=True)
    out.to_csv(dest / "library_composition.tsv", sep="\t", index=False)
    print(out.to_string(index=False))

    if hist_rows:
        h = pd.concat(hist_rows).pivot_table(index="gene", columns="library", values="tpm")
        h.to_csv(dest / "histone_tpm_by_library.tsv", sep="\t")
        print("\nreplication-dependent histones, TPM after rescaling (top 12 by mean):")
        print(h.assign(m=h.mean(axis=1)).nlargest(12, "m").drop(columns="m").round(1).to_string())
    print(f"\nwritten to {dest}")


if __name__ == "__main__":
    main()
