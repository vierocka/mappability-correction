#!/usr/bin/env python3
"""Do demoted (inflated) transcripts carry fewer MS ligands per unit apparent abundance?

The prediction. If pseudogene-derived inflation is an artefact, a transcript whose signal is
inflated has less real protein than its UNCORRECTED abundance implies, so at a given uncorrected
TPM it should yield fewer MS ligands than a transcript whose signal is clean. Unit of analysis is
the transcript, so each observation is independent (peptide-level tests pseudoreplicate: one
protein contributes hundreds of 9-mers).

  outcome    n_ms      MS ligands attributed to the transcript
  exposure   demoted   the correction reduced its counts (f_unique < 1)
  covariates log10 uncorrected TPM, log10 number of 9-mers (protein length)

Assumptions checked before fitting, and why the obvious model is wrong:

  * Ligand yield is NOT proportional to TPM x length. The observed rate per unit exposure falls
    ~145-fold across exposure deciles (detection saturates), so log(TPM x length) must NOT be used
    as an offset; both enter as free covariates, with a quadratic term in log TPM.
  * Counts are overdispersed (var/mean ~ 1.2-1.9) with >90 % zeros, so Poisson standard errors are
    too small. The primary fit is negative binomial (NB2, variance mu + alpha*mu^2) by maximum
    likelihood, with alpha estimated jointly. A Poisson fit with HC0 sandwich errors is reported
    alongside as a sensitivity check: its mean estimate is consistent under any variance
    misspecification, so agreement between the two indicates the result is not an artefact of the
    variance model.
  * Model misspecification is guarded against with a second, model-free test: a Mantel-Haenszel
    rate ratio over (log-TPM ventile x log-length tertile) strata with a within-stratum permutation
    null. It must be a stratum-weighted estimator; a raw difference of means across strata stays
    confounded however finely the strata are cut.

usage: ligand_yield_test.py --data-dir DIR [--outdir .]
"""
import argparse
import pathlib
import warnings

import numpy as np
import pandas as pd
from scipy import optimize, special, stats

from tail_validation import MEAN_COUNTS, MS_FILES, load_ms

warnings.filterwarnings("ignore")
RNG = np.random.default_rng(0)


def build(data, cell):
    t = pd.read_csv(data / cell / "corrected_counts" / MEAN_COUNTS[cell], sep="\t")
    t["f_unique"] = np.where(t.mean_counts > 0, t.mean_counts_corrected / t.mean_counts, np.nan)
    ms = load_ms(data, cell)
    m = pd.read_csv(data / cell / f"master_file_abundance_{cell}.csv",
                    usecols=["transcript_id", "peptide"])
    m["is_ms"] = m.peptide.isin(ms)
    g = m.groupby("transcript_id").agg(n_9mer=("peptide", "size"), n_ms=("is_ms", "sum"))
    d = t.set_index("gene_id").join(g, how="inner")
    d = d[(d.mean_TPM > 0) & (d.n_9mer > 0)].copy()
    d["demoted"] = (d.f_unique < 0.999999).astype(int)
    d["l_tpm"] = np.log10(d.mean_TPM)
    d["l_len"] = np.log10(d.n_9mer)
    d["cell_line"] = cell
    return d


def diagnostics(d, cell):
    v, mu = d.n_ms.var(), d.n_ms.mean()
    e = d.mean_TPM * d.n_9mer
    q = pd.qcut(np.log10(e), 10, labels=False)
    rate = d.groupby(q).apply(lambda x: x.n_ms.sum() / (x.mean_TPM * x.n_9mer).sum(),
                              include_groups=False)
    return {"cell_line": cell, "n_transcripts": len(d), "n_demoted": int(d.demoted.sum()),
            "ligands_total": int(d.n_ms.sum()), "ligands_on_demoted": int(d.loc[d.demoted == 1, "n_ms"].sum()),
            "pct_zero": round(100 * (d.n_ms == 0).mean(), 1),
            "var_over_mean": round(v / mu, 2),
            "rate_decile1_over_decile10": round(float(rate.iloc[0] / rate.iloc[-1]), 1)}


def nb2_mle(y, X):
    """Negative binomial (NB2) by maximum likelihood; returns beta, SE, alpha.

    var = mu + alpha*mu^2. Parameters are (beta, log alpha); SEs come from the numerical
    Hessian of the negative log-likelihood at the optimum.
    """
    X = np.asarray(X, float); y = np.asarray(y, float)

    def nll(p):
        b, la = p[:-1], p[-1]
        al = np.exp(la)
        mu = np.exp(np.clip(X @ b, -30, 30))
        r = 1.0 / al
        return -np.sum(special.gammaln(y + r) - special.gammaln(r) - special.gammaln(y + 1)
                       + r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu)))

    p0 = np.zeros(X.shape[1] + 1)
    p0[0] = np.log(max(y.mean(), 1e-6))
    res = optimize.minimize(nll, p0, method="BFGS",
                            options={"maxiter": 2000, "gtol": 1e-7})
    p = res.x
    # numerical Hessian -> covariance
    n = len(p); H = np.zeros((n, n)); h = 1e-5
    for i in range(n):
        for j in range(i, n):
            pp = p.copy(); pp[i] += h; pp[j] += h; f1 = nll(pp)
            pp = p.copy(); pp[i] += h; pp[j] -= h; f2 = nll(pp)
            pp = p.copy(); pp[i] -= h; pp[j] += h; f3 = nll(pp)
            pp = p.copy(); pp[i] -= h; pp[j] -= h; f4 = nll(pp)
            H[i, j] = H[j, i] = (f1 - f2 - f3 + f4) / (4 * h * h)
    try:
        cov = np.linalg.inv(H)
        se = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        se = np.full(n, np.nan)
    return p[:-1], se[:-1], float(np.exp(p[-1])), bool(res.success)


def poisson_robust(y, X):
    """Poisson GLM by IRLS with HC0 sandwich standard errors (valid under overdispersion)."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    b = np.zeros(X.shape[1])
    for _ in range(200):
        eta = np.clip(X @ b, -30, 30)
        mu = np.exp(eta)
        W = mu
        z = eta + (y - mu) / np.maximum(mu, 1e-12)
        XtW = X.T * W
        b_new = np.linalg.solve(XtW @ X + 1e-10 * np.eye(X.shape[1]), XtW @ z)
        if np.max(np.abs(b_new - b)) < 1e-10:
            b = b_new
            break
        b = b_new
    mu = np.exp(np.clip(X @ b, -30, 30))
    bread = np.linalg.inv((X.T * mu) @ X + 1e-10 * np.eye(X.shape[1]))
    meat = (X.T * (y - mu) ** 2) @ X
    V = bread @ meat @ bread
    se = np.sqrt(np.diag(V))
    disp = float(np.sum((y - mu) ** 2 / np.maximum(mu, 1e-12)) / (len(y) - X.shape[1]))
    return b, se, disp


def fit(d, extra=None):
    cols = {"const": 1.0, "demoted": d.demoted.astype(float), "l_tpm": d.l_tpm,
            "l_tpm2": d.l_tpm ** 2, "l_len": d.l_len}
    if extra:
        cols.update(extra)
    X = pd.DataFrame(cols)
    i = list(X.columns).index("demoted")
    y = d.n_ms.values

    b, se, alpha, ok = nb2_mle(y, X.values)
    z = b[i] / se[i]
    out = {"alpha": round(alpha, 3), "converged": ok,
           "RR_nb": round(float(np.exp(b[i])), 3),
           "ci_low_nb": round(float(np.exp(b[i] - 1.96 * se[i])), 3),
           "ci_high_nb": round(float(np.exp(b[i] + 1.96 * se[i])), 3),
           "p_nb": float(2 * stats.norm.sf(abs(z)))}

    bp, sep, disp = poisson_robust(y, X.values)
    zp = bp[i] / sep[i]
    out.update({"dispersion": round(disp, 2),
                "RR_poisson_robust": round(float(np.exp(bp[i])), 3),
                "p_poisson_robust": float(2 * stats.norm.sf(abs(zp)))})
    return out


def matched(d, n_perm=20000, n_tpm=20):
    """Model-free stratified test: Mantel-Haenszel rate ratio with a within-stratum permutation null.

    Strata are log-TPM ventiles x log-length tertiles. Within each stratum the ligands expected on
    demoted transcripts under the null are E_s = n_demoted_s * (ligands_s / n_s); the estimator is
    sum(O_s) / sum(E_s), which is stratum-weighted and so directly comparable to the GLM rate ratio.
    A raw difference of means across strata is NOT comparable: it is a marginal contrast and stays
    confounded however finely the strata are cut.
    """
    d = d.copy()
    d["s"] = (pd.qcut(d.l_tpm, n_tpm, labels=False, duplicates="drop").astype(str) + "_" +
              pd.qcut(d.l_len, 3, labels=False, duplicates="drop").astype(str))
    keep = d.groupby("s").demoted.transform(lambda x: 0 < x.sum() < len(x))
    d = d[keep]
    if d.empty or d.demoted.sum() == 0:
        return {"n_strata": 0}

    y = d.n_ms.values.astype(float)
    lab = d.demoted.values.astype(int)
    strata = d.s.values
    idx = [np.flatnonzero(strata == s) for s in np.unique(strata)]

    def mh(labels):
        O = E = 0.0
        for ii in idx:
            k = labels[ii].sum()
            if k:
                O += y[ii][labels[ii] == 1].sum()
                E += k * y[ii].mean()
        return O, E

    O, E = mh(lab)
    rr = O / E if E > 0 else np.nan
    rng = np.random.default_rng(0)
    null = np.empty(n_perm)
    for b in range(n_perm):
        p = lab.copy()
        for ii in idx:
            p[ii] = rng.permutation(lab[ii])
        o, e = mh(p)
        null[b] = o / e if e > 0 else np.nan
    pv = (np.sum(np.abs(np.log(null)) >= abs(np.log(rr))) + 1) / (n_perm + 1)
    lo, hi = np.exp(np.log(rr) - 1.96 * np.std(np.log(null))), np.exp(np.log(rr) + 1.96 * np.std(np.log(null)))
    return {"n_strata": len(idx), "n_demoted_matched": int(lab.sum()),
            "ligands_obs": round(float(O), 1), "ligands_exp": round(float(E), 1),
            "RR_mh": round(float(rr), 3), "mh_ci_low": round(float(lo), 3),
            "mh_ci_high": round(float(hi), 3), "perm_p": round(float(pv), 4)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data-dir", required=True)
    p.add_argument("--cells", nargs="+", default=["JY", "HCC1937", "HCT116"])
    p.add_argument("--outdir", default=".")
    a = p.parse_args()
    data = pathlib.Path(a.data_dir).expanduser()
    out = pathlib.Path(a.outdir)

    diag, res, frames = [], [], []
    for cell in a.cells:
        d = build(data, cell)
        frames.append(d)
        diag.append(diagnostics(d, cell))
        r = {"cell_line": cell}
        r.update(fit(d))
        r.update(matched(d))
        res.append(r)

    pooled = pd.concat(frames)
    rp = {"cell_line": "POOLED"}
    rp.update(fit(pooled, {f"cell_{c}": (pooled.cell_line == c).astype(float)
                           for c in a.cells[1:]}))
    rp.update(matched(pooled))
    res.append(rp)

    dd, rr = pd.DataFrame(diag), pd.DataFrame(res)
    dd.to_csv(out / "ligand_yield_diagnostics.tsv", sep="\t", index=False)
    rr.to_csv(out / "ligand_yield_test.tsv", sep="\t", index=False)
    print("ASSUMPTION DIAGNOSTICS")
    print(dd.to_string(index=False))
    print("\nNEGATIVE BINOMIAL (NB2, ML) — rate ratio for 'demoted', adjusted for log TPM, its square, log length")
    print(rr[["cell_line", "alpha", "RR_nb", "ci_low_nb", "ci_high_nb", "p_nb", "dispersion", "RR_poisson_robust", "p_poisson_robust"]].to_string(index=False))
    print("\nMODEL-FREE STRATIFIED (Mantel-Haenszel RR, within-stratum permutation null)")
    print(rr[["cell_line", "n_strata", "n_demoted_matched", "ligands_obs", "ligands_exp",
              "RR_mh", "mh_ci_low", "mh_ci_high", "perm_p"]].to_string(index=False))
    print("\nrate ratio < 1 = demoted transcripts carry FEWER ligands than their apparent abundance implies")


if __name__ == "__main__":
    main()
