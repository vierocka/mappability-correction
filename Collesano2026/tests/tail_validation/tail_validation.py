#!/usr/bin/env python3
"""Does the mappability correction improve agreement with MS *where it acts*?

Global recall and AUC cannot answer this: the correction reaches ~0.1 % of the 9.4 M 9-mers, so a
rank statistic over all of them is dominated by the 99.9 % it never touches. This test restricts
the comparison to the expression tail, which is where the correction is concentrated
(31–36 % of the top 0.1 % of net sums change; see ../inflation_tail/).

Net sums are rebuilt twice from the same peptide/transcript table:
    before = sum over occurrences of the UNCORRECTED transcript TPM
    after  = sum over occurrences of the CORRECTED   transcript TPM
so the two rankings differ only by x funique and the exclusions.

Reported at each tail depth k:
    PPV@k     of the top k 9-mers by net sum, the fraction that are MS-detected ligands
    recall@k  of the MS ligands, the fraction that fall in the top k
A correction that sharpens the tail raises both. The absolute values are low by construction:
the source repertoire has no HLA-binding step, so most high-abundance 9-mers are not ligands.
Only the before/after difference is interpretable.

Inputs are read-only; nothing is written outside this folder.

usage: tail_validation.py --data-dir DIR [--cells JY HCC1937 HCT116] [--outdir .]
"""
import argparse
import pathlib

import numpy as np
import pandas as pd

MEAN_COUNTS = {  # per cell line: the replicate-mean corrected-counts table
    "JY": "JY_corrected_mean_364065_mane_fc_corrected.tsv",
    "HCC1937": "HCC1937_corrected_mean_74_75_mane_fc_corrected.tsv",
    "HCT116": "HCT116_corrected_mean_41_42_mane_fc_corrected.tsv",
}
MS_FILES = {  # MS-detected ligands; JY is typed for two alleles and both are used
    "JY": ["MS_total_expression_JY0201.csv", "MS_total_expression_JY0702.csv"],
    "HCC1937": ["MS_total_expression_HCC1937.csv"],
    "HCT116": ["MS_total_expression_HCT116.csv"],
}
DEPTHS = [0.001, 0.005, 0.01, 0.05, 0.10]   # top 0.1 %, 0.5 %, 1 %, 5 %, 10 %


def load_ms(data, cell):
    pep = set()
    for fn in MS_FILES[cell]:
        d = pd.read_csv(data / cell / fn, usecols=["9mer"])
        pep |= set(d["9mer"].astype(str))
    return pep


def net_sums(data, cell):
    """Per-9-mer net sums before and after correction, plus the number of source transcripts."""
    t = pd.read_csv(data / cell / "corrected_counts" / MEAN_COUNTS[cell], sep="\t",
                    usecols=["gene_id", "mean_TPM", "mean_TPM_corrected"]).set_index("gene_id")
    m = pd.read_csv(data / cell / f"master_file_abundance_{cell}.csv",
                    usecols=["transcript_id", "peptide"])
    m["before"] = m.transcript_id.map(t.mean_TPM)
    m["after"] = m.transcript_id.map(t.mean_TPM_corrected)
    n_unmapped = int(m.before.isna().sum())
    m = m[m.before.notna()]
    g = m.groupby("peptide")
    pep = pd.DataFrame({"before": g.before.sum(), "after": g.after.sum(),
                        "n_src": g.transcript_id.nunique()})
    pep = pep[pep.before > 0]
    return pep, n_unmapped


def tail_stats(pep, ms, col, k):
    n = max(int(round(len(pep) * k)), 1)
    top = pep.nlargest(n, col)
    hits = top.index.isin(ms).sum()
    return n, hits, 100 * hits / n, 100 * hits / len(ms)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data-dir", required=True)
    p.add_argument("--cells", nargs="+", default=["JY", "HCC1937", "HCT116"])
    p.add_argument("--outdir", default=".")
    a = p.parse_args()
    data = pathlib.Path(a.data_dir).expanduser()
    out = pathlib.Path(a.outdir)

    rows, rank_rows = [], []
    for cell in a.cells:
        ms = load_ms(data, cell)
        pep, n_unmapped = net_sums(data, cell)
        in_rep = len(ms & set(pep.index))
        print(f"{cell}: {len(pep):,} 9-mers, {len(ms):,} MS ligands, "
              f"{in_rep} of them in the repertoire ({100*in_rep/len(ms):.1f} % global recall); "
              f"{n_unmapped:,} occurrences dropped (transcript absent from the corrected table)")
        # every TPM shifts slightly because correcting 462 transcripts renormalises the library;
        # the meaningful count is how many move RELATIVE to that common rescaling
        shift = np.log2(pep.after / pep.before)
        base = shift.median()
        pep["rel_shift"] = shift - base
        moved = int((pep.rel_shift.abs() > 1e-6).sum())
        print(f"   renormalisation factor {2**base:.4f}x; {moved:,} 9-mers ({100*moved/len(pep):.2f} %) "
              f"move relative to it")
        for k in DEPTHS:
            r = {"cell_line": cell, "depth_pct": 100 * k, "n_ms": len(ms)}
            for col in ("before", "after"):
                n, hits, ppv, rec = tail_stats(pep, ms, col, k)
                r[f"n_top"] = n
                r[f"hits_{col}"] = hits
                r[f"ppv_{col}"] = round(ppv, 4)
                r[f"recall_{col}"] = round(rec, 3)
            r["d_ppv"] = round(r["ppv_after"] - r["ppv_before"], 4)
            r["d_recall"] = round(r["recall_after"] - r["recall_before"], 3)
            rows.append(r)

        # paired rank test: 19 MS hits in the top 0.1 % is too few to compare PPV on, but all
        # MS ligands can be followed individually. Rank 1 = highest net sum.
        rb = pep.before.rank(ascending=False, method="average")
        ra = pep.after.rank(ascending=False, method="average")
        idx = pep.index.isin(ms)
        db, da = rb[idx], ra[idx]
        delta = (db - da)                      # positive = ligand moved UP after correction
        up, down = int((delta > 0).sum()), int((delta < 0).sum())
        # Matched control. Correcting some peptides downward lifts everything else in rank, so
        # "ligands moved up" is expected under the null. The question is whether they move up
        # MORE than non-ligands of the same abundance, which is what the correction would predict.
        bins = pd.qcut(np.log10(pep.before.clip(lower=1e-9)), 100, labels=False, duplicates="drop")
        rng = np.random.default_rng(0)
        ctrl = []
        for b, n_need in pd.Series(bins[idx]).value_counts().items():
            pool = np.flatnonzero((bins.values == b) & ~idx)
            if len(pool):
                ctrl.append(rng.choice(pool, size=min(n_need * 10, len(pool)), replace=False))
        ctrl = np.concatenate(ctrl) if ctrl else np.array([], dtype=int)
        d_ctrl = (rb.values - ra.values)[ctrl]
        med_ms, med_ctrl = float(np.median(delta)), float(np.median(d_ctrl))
        try:
            from scipy.stats import mannwhitneyu
            pv_ctrl = mannwhitneyu(delta, d_ctrl, alternative="two-sided").pvalue
        except Exception:
            pv_ctrl = float("nan")
        print(f"   vs abundance-matched non-ligands (n={len(ctrl):,}): "
              f"ligand median {med_ms:+.0f} vs control {med_ctrl:+.0f}  (Mann-Whitney p = {pv_ctrl:.3g})")
        try:
            from scipy.stats import wilcoxon
            pv = wilcoxon(db, da, zero_method="zsplit").pvalue if up + down else float("nan")
        except Exception:
            pv = float("nan")
        print(f"   MS ligands: {up} move up, {down} move down, median rank change "
              f"{np.median(delta):+.0f} (Wilcoxon p = {pv:.3g})")
        rank_rows.append({"cell_line": cell, "n_ms": int(idx.sum()), "moved_up": up,
                          "moved_down": down, "median_rank_change": float(np.median(delta)),
                          "median_rel_shift_ms": float(pep.rel_shift[idx].median()),
                          "median_rel_shift_all": float(pep.rel_shift.median()),
                          "wilcoxon_p_no_null": pv,
                          "median_rank_change_control": med_ctrl,
                          "n_control": int(len(ctrl)),
                          "mannwhitney_p_vs_control": pv_ctrl})

    d = pd.DataFrame(rows)
    d.to_csv(out / "tail_validation.tsv", sep="\t", index=False)
    pd.DataFrame(rank_rows).to_csv(out / "ms_rank_shift.tsv", sep="\t", index=False)
    print()
    print(d.to_string(index=False))
    print()
    print("positive d_ppv / d_recall = the correction improves agreement with MS in that tail")


if __name__ == "__main__":
    main()
