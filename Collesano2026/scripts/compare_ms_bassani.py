#!/usr/bin/env python3
"""
Predicted 9-mer peptidome vs HLA-I ligands of Bassani-Sternberg 2015 (MCP, PMID 25576301).
MS truth = Table S1 (MaxQuant Peptides), Intensity > 0 in any replicate of the cell line.

Per MS 9-mer, one fate:
  caught            in the source repertoire (net sum TPM_corrected > 0)
  not_expressed     in MANE proteome, source transcript(s) TPM_corrected = 0
  filtered_out      in MANE proteome, source transcript(s) excluded (NaN)
  not_in_MANE       absent from MANE proteome; split by presence in the study's UniProt DB
Metrics: recall of MS ligands, ROC-AUC of raw vs corrected net-sum ranks.
The SOURCE REPERTOIRE is the set of 9-mers enumerated from the expressed personal MANE
proteome (net sum > 0). It is not a binding prediction: no HLA binding, processing or
length model is applied. Jaccard against the ligandome is not reported, being a ratio of
set sizes (the source repertoire exceeds the ligandome ~3000-fold), not a measure of agreement.
"""
import argparse
import re

import numpy as np
import pandas as pd


def auc(score, label):
    """Mann-Whitney AUC: P(score_MS > score_nonMS)."""
    r = pd.Series(score).rank().values
    n1 = label.sum(); n0 = len(label) - n1
    return (r[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def uniprot_9mers(fasta, wanted):
    seqs, cur = [], []
    for line in open(fasta):
        if line.startswith(">"):
            if cur: seqs.append("".join(cur))
            cur = []
        else:
            cur.append(line.strip())
    if cur: seqs.append("".join(cur))
    hit = set()
    for s in seqs:
        for p in wanted - hit:
            if p in s: hit.add(p)
    return hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table-s1", required=True)
    ap.add_argument("--uniprot-fasta", required=True, help="study search DB (HUMAN.fasta)")
    ap.add_argument("--cell", required=True, help="column stem in Table S1, e.g. JY")
    ap.add_argument("--peptides", required=True, help="peptide_abundance_<CELL>.csv")
    ap.add_argument("--master", required=True, help="master_file_abundance_<CELL>.csv")
    ap.add_argument("--k", type=int, default=9)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    s1 = pd.read_excel(a.table_s1, sheet_name="Peptides")
    cols = [c for c in s1.columns if re.fullmatch(rf"Intensity {re.escape(a.cell)}-\d+", c)]
    ms_all = s1.loc[s1[cols].fillna(0).gt(0).any(axis=1), ["Sequence", "Length", "Gene names"]]
    len_dist = ms_all.Length.value_counts().sort_index()
    ms = set(ms_all.loc[ms_all.Length == a.k, "Sequence"])

    pep = pd.read_csv(a.peptides, usecols=["peptide", "mean_TPM", "mean_TPM_corrected"])
    source_repertoire = set(pep.loc[pep.mean_TPM_corrected > 0, "peptide"])

    # MANE proteome membership + source status (all rows, incl. excluded / zero)
    src = []
    for ch in pd.read_csv(a.master, usecols=["transcript_id", "gene_symbol", "peptide", "mean_TPM_corrected"],
                          chunksize=2_000_000):
        src.append(ch[ch.peptide.isin(ms)])
    src = pd.concat(src)
    best = src.groupby("peptide").mean_TPM_corrected.agg(lambda x: x.max(skipna=True) if x.notna().any() else np.nan)

    fate = {}
    for p in ms:
        if p in source_repertoire:
            fate[p] = "caught"
        elif p in best.index:
            fate[p] = "filtered_out" if np.isnan(best[p]) else "not_expressed"
        else:
            fate[p] = "not_in_MANE"
    miss = {p for p, f in fate.items() if f == "not_in_MANE"}
    in_uni = uniprot_9mers(a.uniprot_fasta, miss)
    for p in miss:
        fate[p] = "not_in_MANE: other UniProt isoform/protein" if p in in_uni else "not_in_MANE: not in UniProt DB"

    f = pd.Series(fate, name="fate").rename_axis("peptide").reset_index()
    f = f.merge(src.groupby("peptide").agg(genes=("gene_symbol", lambda x: ",".join(sorted(set(x))))),
                on="peptide", how="left")
    f = f.merge(pep, on="peptide", how="left")
    f.to_csv(f"{a.out_prefix}_ms_peptide_fates.tsv", sep="\t", index=False)

    inter = len(ms & source_repertoire)
    pp = pep[pep.mean_TPM_corrected > 0].copy()
    lab = pp.peptide.isin(ms).values
    summ = {
        "cell_line": a.cell, "ms_peptides_all_lengths": len(ms_all),
        f"ms_{a.k}mers": len(ms), f"pct_ms_{a.k}mers": round(100 * len(ms) / len(ms_all), 1),
        "source_repertoire_9mers": len(source_repertoire), "overlap": inter,
        "recall_pct": round(100 * inter / len(ms), 2),
        "auc_raw_TPM": round(auc(pp.mean_TPM.values, lab), 4),
        "auc_corrected_TPM": round(auc(pp.mean_TPM_corrected.values, lab), 4),
    }
    for k, v in f.fate.value_counts().items():
        summ[f"n_{k}"] = v
        summ[f"pct_{k}"] = round(100 * v / len(ms), 2)
    pd.Series(summ).to_frame("value").to_csv(f"{a.out_prefix}_summary.tsv", sep="\t")
    len_dist.rename("n").to_csv(f"{a.out_prefix}_ms_length_distribution.tsv", sep="\t")

    # MS share per expression decile (corrected TPM)
    pp["decile"] = pd.qcut(pp.mean_TPM_corrected.rank(method="first"), 10, labels=range(1, 11))
    pp["is_ms"] = lab
    dec = pp.groupby("decile", observed=True).agg(n=("peptide", "size"), n_ms=("is_ms", "sum"),
                                                   tpm_min=("mean_TPM_corrected", "min"))
    dec["ms_per_10k"] = 1e4 * dec.n_ms / dec.n
    dec.to_csv(f"{a.out_prefix}_ms_by_expression_decile.tsv", sep="\t")
    print(pd.Series(summ).to_string())


if __name__ == "__main__":
    main()
