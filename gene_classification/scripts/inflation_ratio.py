#!/usr/bin/env python3
"""Inflation: what a UTR-containing reference attributes that a CDS-only reference cannot.

The same simulated reads are pseudoaligned against two references: one transcript per gene with
UTRs, and the deduplicated CDS set. Reads are uniform across every transcript by construction,
so any difference between the two TPM estimates is attribution, not expression.

inflation_ratio = TPM(representative, UTRs) / TPM(CDS, deduplicated), per gene.

A high ratio means UTR k-mers let the pseudoaligner place reads that the coding sequence alone
cannot resolve. That is useful when the UTR is genuinely unique, and misleading when the UTR is
shared while the CDS is not - so the ratio is reported, not acted on by itself.

usage: inflation_ratio.py --utr F --cds F --tx2gene F --out F [--pseudo 0.01]
"""
import argparse

import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--utr", required=True)
    p.add_argument("--cds", required=True)
    p.add_argument("--tx2gene", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--pseudo", type=float, default=0.01, help="floor for the ratio denominator")
    p.add_argument("--inflation-hi", type=float, default=2.0,
                   help="ratio at or above which a gene counts as UTR-inflated")
    a = p.parse_args()

    t2g = pd.read_csv(a.tx2gene, sep="\t")
    name = dict(zip(t2g.gene_id, t2g.gene_name))

    utr = pd.read_csv(a.utr, sep="\t")[["target_id", "tpm"]]
    utr["transcript_id"] = utr.target_id.astype(str).str.split().str[0]
    utr = utr.merge(t2g[["transcript_id", "gene_id"]], on="transcript_id", how="left")
    utr_g = utr.groupby("gene_id", as_index=False).tpm.sum().rename(columns={"tpm": "tpm_utr"})

    # CDS entries are not transcript-named in every annotation; map by gene where possible,
    # otherwise by the transcript id embedded in the CDS name.
    cds = pd.read_csv(a.cds, sep="\t")[["target_id", "tpm"]]
    cds["transcript_id"] = cds.target_id.astype(str).str.split().str[0]
    cds = cds.merge(t2g[["transcript_id", "gene_id"]], on="transcript_id", how="left")
    cds_g = cds.groupby("gene_id", as_index=False).tpm.sum().rename(columns={"tpm": "tpm_cds"})

    d = utr_g.merge(cds_g, on="gene_id", how="outer").fillna({"tpm_utr": 0.0, "tpm_cds": 0.0})
    d["gene_name"] = d.gene_id.map(name)
    d["inflation_ratio"] = d.tpm_utr / d.tpm_cds.clip(lower=a.pseudo)
    d.loc[(d.tpm_utr == 0) & (d.tpm_cds == 0), "inflation_ratio"] = pd.NA

    # plain strings, not pd.cut: a Categorical rejects the extra no_cds_signal label
    def label(r):
        if r.tpm_cds <= a.pseudo:
            return "no_cds_signal"
        if pd.isna(r.inflation_ratio):
            return "no_signal"
        if r.inflation_ratio < 0.5:
            return "cds_inflated"
        if r.inflation_ratio < a.inflation_hi:
            return "consistent"
        return "utr_inflated"

    d["class"] = d.apply(label, axis=1)

    d[["gene_id", "gene_name", "tpm_utr", "tpm_cds", "inflation_ratio", "class"]] \
        .to_csv(a.out, sep="\t", index=False)

    print(f"{len(d)} genes")
    print(d["class"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
