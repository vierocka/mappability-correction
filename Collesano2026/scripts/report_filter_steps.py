#!/usr/bin/env python3
"""
Per sample: what each filtering / correction step removes, and how the library
composition shifts when excluded transcripts leave the TPM denominator.
Input: <SRR>_mane_fc_corrected.tsv (apply_uniqueness_and_average.py output).

Outputs (prefix):
  _steps.tsv        step, n_transcripts, counts_removed, % raw counts, % raw TPM, % corrected TPM
  _groups.tsv       transcripts per correction rule (G0-G7 of README_decisions.md)
  _composition.tsv  raw TPM share by class (non-coding by biotype, PAR, alt, MT, retained)
  _shift.tsv        log2(TPM_corrected / TPM) of uncorrected transcripts: median, IQR
  _top_excluded.tsv 15 excluded transcripts with the largest raw TPM share
  _*.png            plots
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def classify(rule, tid, fix_ids):
    if rule.startswith("excluded: non-coding"): return "non_coding"
    if rule.startswith("excluded: PAR"): return "PAR"
    if rule.startswith("excluded: alt copy"): return "alt_copy_folded"
    if rule.startswith("excluded: non-primary"):
        return "no_unique_locus" if tid in fix_ids else "alt_population_variable"
    return "retained"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True, help="<SRR>_mane_fc_corrected.tsv files")
    ap.add_argument("--labels", nargs="+", required=True, help="cell_line:sample per file")
    ap.add_argument("--resources", required=True)
    ap.add_argument("--min-tpm", type=float, default=1.0)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    R = a.resources.rstrip("/")
    npm = pd.read_csv(f"{R}/non_primary_MANE_transcripts.tsv", sep="\t")
    fix_ids = set(npm.loc[npm.gff_chromosome.str.endswith("_fix"), "transcript_id"])
    nc = pd.read_csv(f"{R}/non_coding_MANE_transcripts_no_protein.tsv", sep="\t")
    biotype = dict(nc[["transcript_id", "transcript_biotype"]].values)
    symbol = dict(nc[["transcript_id", "gene_name"]].values) | dict(npm[["transcript_id", "gene_name"]].values)

    steps, groups, comp, shift, top = [], [], [], [], []
    for path, lab in zip(a.samples, a.labels):
        d = pd.read_csv(path, sep="\t")
        d["rule"] = d.rule.fillna("")
        d["cls"] = [classify(r, t.split("_")[0], fix_ids) for r, t in zip(d.rule, d.gene_id)]
        C, T = d.counts.sum(), d.TPM.sum()
        for cls in ["non_coding", "PAR", "alt_copy_folded", "alt_population_variable", "no_unique_locus"]:
            x = d[d.cls == cls]
            steps.append([lab, f"excluded: {cls}", len(x), x.counts.sum(), 100 * x.counts.sum() / C,
                          100 * x.TPM.sum() / T, np.nan,
                          "counts moved to primary copy" if cls == "alt_copy_folded" else "removed"])
        r3 = d[d.rule.str.startswith("pseudogene inflation")]
        lost = (r3.counts - r3.counts_corrected).sum()
        steps.append([lab, "x funique (pseudogene inflation)", len(r3), lost, 100 * lost / C, np.nan, np.nan,
                      f"median f = {r3.f_unique.median():.3f}"])
        kept = d[d.cls == "retained"]
        low = kept[kept.TPM_corrected < a.min_tpm]
        steps.append([lab, f"TPM_corrected < {a.min_tpm:g}", len(low), low.counts_corrected.sum(),
                      100 * low.counts_corrected.sum() / C, 100 * low.TPM.sum() / T,
                      100 * low.TPM_corrected.sum() / 1e6, "removed before peptides"])
        fin = kept[kept.TPM_corrected >= a.min_tpm]
        steps.append([lab, "final (peptide input)", len(fin), fin.counts_corrected.sum(),
                      100 * fin.counts_corrected.sum() / C, 100 * fin.TPM.sum() / T,
                      100 * fin.TPM_corrected.sum() / 1e6, "kept"])

        g = kept.rule.str.replace(r";.*", "", regex=True).str.replace(r" \(\+alt fcM\)", "", regex=True)
        for k, v in g.value_counts().items():
            groups.append([lab, k, v])

        cl = np.where(d.cls == "non_coding", "non_coding: " + d.gene_id.map(biotype).fillna("?"), d.cls)
        s = d.assign(cls2=cl).groupby("cls2").agg(raw_TPM_pct=("TPM", lambda x: 100 * x.sum() / T),
                                                   raw_counts_pct=("counts", lambda x: 100 * x.sum() / C),
                                                   n=("gene_id", "size"))
        for k, r in s.iterrows():
            comp.append([lab, k, r.n, r.raw_counts_pct, r.raw_TPM_pct])

        u = kept[(kept.rule.str.startswith("f=1")) & (kept.TPM > 0) & (kept.TPM_corrected > 0)]
        l2 = np.log2(u.TPM_corrected / u.TPM)
        shift.append([lab, len(u), l2.median(), l2.quantile(.25), l2.quantile(.75)])
        ex = d[d.cls != "retained"].nlargest(15, "TPM")
        for _, r in ex.iterrows():
            top.append([lab, r.gene_id, symbol.get(r.gene_id, ""), r.cls, 100 * r.TPM / T])

    P = a.out_prefix
    st = pd.DataFrame(steps, columns=["sample", "step", "n_transcripts", "counts", "pct_raw_counts",
                                      "pct_raw_TPM", "pct_corrected_TPM", "note"])
    st.to_csv(f"{P}_steps.tsv", sep="\t", index=False, float_format="%.3f")
    pd.DataFrame(groups, columns=["sample", "rule", "n_transcripts"]).to_csv(f"{P}_groups.tsv", sep="\t", index=False)
    co = pd.DataFrame(comp, columns=["sample", "class", "n_transcripts", "pct_raw_counts", "pct_raw_TPM"])
    co.to_csv(f"{P}_composition.tsv", sep="\t", index=False, float_format="%.3f")
    sh = pd.DataFrame(shift, columns=["sample", "n_uncorrected", "median_log2_shift", "q25", "q75"])
    sh.to_csv(f"{P}_shift.tsv", sep="\t", index=False, float_format="%.4f")
    pd.DataFrame(top, columns=["sample", "transcript_id", "gene", "class", "pct_raw_TPM"]).to_csv(
        f"{P}_top_excluded.tsv", sep="\t", index=False, float_format="%.3f")

    # plot 1: raw TPM composition (excluded classes only; retained = remainder)
    pv = co[co["class"] != "retained"].pivot(index="sample", columns="class", values="pct_raw_TPM").fillna(0)
    ax = pv.plot.bar(stacked=True, figsize=(7, 4), width=.7, colormap="tab20")
    ax.set_ylabel("% of raw TPM (excluded classes)"); ax.set_xlabel("")
    ax.set_title("Library share removed by exclusion"); ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1, 1))
    plt.xticks(rotation=30, ha="right"); plt.tight_layout(); plt.savefig(f"{P}_composition.png", dpi=150); plt.close()

    # plot 2: step losses (% raw counts)
    pv = st[~st.step.str.startswith("final")].pivot(index="step", columns="sample", values="pct_raw_counts")
    ax = pv.plot.barh(figsize=(7, 4)); ax.set_xlabel("% of raw counts"); ax.set_ylabel("")
    ax.set_title("Counts removed per step"); ax.legend(fontsize=7, frameon=False)
    plt.tight_layout(); plt.savefig(f"{P}_steps.png", dpi=150); plt.close()
    print(st.to_string(index=False)); print(sh.to_string(index=False))


if __name__ == "__main__":
    main()
