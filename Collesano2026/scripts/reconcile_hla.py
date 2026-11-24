#!/usr/bin/env python3
"""
Reconcile OptiType vs arcasHLA class-I calls (2-field) and check against Cellosaurus.

Per gene (A, B, C):
  1. Candidates = union of OptiType, arcasHLA and Cellosaurus alleles (IMGT CDS, select_hla_alleles rule).
  2. HLA reads (fished) are aligned to each candidate CDS separately (bwa mem); per read the best AS.
  3. Genotype support = reads explained by the pair (best AS of the pair == best AS over all candidates
     of that gene); penalty = sum of AS lost by the pair.
  4. Site check: at MAFFT columns where disputed alleles differ, an allele whose own base is seen in
     < 15 % (median) or < 5 % at any site with >= 10 reads is rejected; only reads whose best hit
     is that gene are counted (e.g. JY B*07:412: 7/204 reads, mean Q 14).
  5. Decision: --override (user-typed) > tools agree > pair with higher support; near tie
     (< max(5 reads, 1 %)) -> Cellosaurus best record > OptiType, flagged REVIEW.
     Cellosaurus: every record scored by read support; best-supported record
     reported as concordance (conflicting/dubious records visible), never silently applied.
  6. Disagreeing genes: MAFFT of candidate CDS -> variable sites -> per site, reads carrying each
     allele's base and their mean base quality; plot.
Outputs (--outdir): <s>_hla_calls.tsv (final), <s>_hla_support.tsv, <s>_hla_sites.tsv,
                    <s>_hla_cellosaurus.tsv, <s>_hla_disagreement.png
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import urllib.parse
import urllib.request
from collections import defaultdict
from itertools import combinations_with_replacement
from pathlib import Path

import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).parent))
from select_hla_alleles import read_nuc, pick, CDS_LEN  # noqa: E402

GENES = ["A", "B", "C"]
MIN_SITE_FRAC = 0.15     # allele-specific base in < 15 % of covering reads -> allele rejected
MIN_SITE_FRAC_ANY = 0.05 # ... or < 5 % at any specific site with >= 10 reads
MIN_SITE_COV = 10
MIN_DELTA_READS = 5      # support difference below max(5 reads, 1 % gene reads) = near tie
MIN_DELTA_PCT = 1.0
API = "https://api.cellosaurus.org"


def two(a):
    g, f = a.replace("HLA-", "").split("*")
    return f"{g}*{':'.join(f.split(':')[:2])}"


def optitype_calls(path):
    t = pd.read_csv(path, sep="\t", index_col=0).iloc[0]
    return {g: sorted([two(t[f"{g}1"]), two(t[f"{g}2"])]) for g in GENES}, float(t["Reads"])


def arcas_calls(path):
    j = json.load(open(path))
    return {g: sorted(two(a) for a in (j[g] * 2)[:2]) for g in GENES if g in j}


# ── Cellosaurus ──────────────────────────────────────────────────────────────
def cellosaurus_lines(query):
    ac = query if re.fullmatch(r"CVCL_\w{4}", query) else None
    if ac is None:
        q = urllib.parse.quote(f'id:"{query}"')
        js = json.load(urllib.request.urlopen(f"{API}/search/cell-line?q={q}&fields=ac,id&format=json", timeout=30))
        hits = [c for c in js["Cellosaurus"]["cell-line-list"]
                if any(n["value"] == query for n in c["name-list"])]
        if not hits:
            raise SystemExit(f"Cellosaurus: no exact identifier '{query}'; pass CVCL_ accession")
        ac = hits[0]["accession-list"][0]["value"]
    txt = urllib.request.urlopen(f"{API}/cell-line/{ac}?format=txt", timeout=30).read().decode()
    return ac, [l[len("CC   HLA typing: "):] for l in txt.splitlines() if l.startswith("CC   HLA typing:")]


def parse_cellosaurus(line):
    """'A*01:01,02:01; B*45:01; ... (PubMed=..)' -> {gene: [a1,a2]} (2-field; 1 allele = hom.)."""
    src = re.search(r"\(([^)]*)\)\s*\.?$", line)
    out = {}
    for part in re.split(r";\s*", re.sub(r"\([^)]*\)\s*\.?$", "", line)):
        m = re.match(r"\s*([ABC])\*(.+)", part)
        if not m:
            continue
        al = [f"{m.group(1)}*{':'.join(x.strip().split(':')[:2])}" for x in m.group(2).split(",") if x.strip()]
        out[m.group(1)] = sorted((al * 2)[:2])
    return out, (src.group(1) if src else "")


def cellosaurus_records(args):
    if args.cellosaurus_file:                              # user-typed: gene allele1 allele2 source
        t = pd.read_csv(args.cellosaurus_file, sep="\t")
        return "user", [({r.gene: sorted([two(r.allele1), two(r.allele2)])}, r.get("source", "user"))
                        for _, r in t.iterrows()]
    if not args.cellosaurus:
        return None, []
    try:
        ac, lines = cellosaurus_lines(args.cellosaurus)
    except Exception as e:                                 # offline: report, continue
        print(f"Cellosaurus lookup failed ({e}); continuing without", file=sys.stderr)
        return None, []
    return ac, [parse_cellosaurus(l) for l in lines]


# ── read support ─────────────────────────────────────────────────────────────
def write_fa(path, name, seq):
    open(path, "w").write(f">{name}\n{seq}\n")


def align_each(cands, reads, tmp, threads):
    """bwa mem reads -> each candidate separately. Returns {allele: bam path}."""
    bams = {}
    for al, (_, seq) in cands.items():
        tag = re.sub(r"\W", "_", al)
        fa = f"{tmp}/{tag}.fa"
        write_fa(fa, al, seq)
        subprocess.run(["bwa", "index", fa], check=True, capture_output=True)
        bam = f"{tmp}/{tag}.bam"
        p = subprocess.Popen(["bwa", "mem", "-t", str(threads), fa, reads], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL)
        subprocess.run(["samtools", "sort", "-o", bam, "-"], stdin=p.stdout, check=True)
        subprocess.run(["samtools", "index", bam], check=True)
        bams[al] = bam
    return bams


def best_as(bams):
    """{read: {allele: AS}} from primary alignments."""
    s = defaultdict(dict)
    for al, bam in bams.items():
        for r in pysam.AlignmentFile(bam):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            s[r.query_name][al] = max(r.get_tag("AS"), s[r.query_name].get(al, -1))
    return s


def genotype_support(scores, cand_by_gene):
    rows, per_gene_reads = [], defaultdict(dict)
    for rd, d in scores.items():                            # assign read to the gene of its best hit
        top = max(d.values())
        genes = {a.split("*")[0] for a, v in d.items() if v == top}
        if len(genes) == 1:
            per_gene_reads[genes.pop()][rd] = d
    for g, al in cand_by_gene.items():
        reads = per_gene_reads[g]
        for x, y in combinations_with_replacement(sorted(al), 2):
            expl = lost = 0
            for d in reads.values():
                top = max(d.values())
                v = max(d.get(x, 0), d.get(y, 0))
                expl += v == top
                lost += top - v
            rows.append({"gene": g, "allele1": x, "allele2": y, "gene_reads": len(reads),
                         "explained": expl, "pct_explained": round(100 * expl / max(len(reads), 1), 2),
                         "AS_lost": lost})
    return pd.DataFrame(rows), per_gene_reads


def mafft(seqs, tmp, g):
    fa = f"{tmp}/{g}_cands.fa"
    open(fa, "w").write("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    out = subprocess.run(["mafft", "--quiet", "--auto", fa], capture_output=True, text=True, check=True).stdout
    aln, h = {}, None
    for ln in out.splitlines():
        if ln.startswith(">"):
            h = ln[1:].strip(); aln[h] = ""
        else:
            aln[h] += ln.strip().upper()
    return aln


def site_evidence(aln, bams, keep=None):
    """Variable MAFFT columns -> per allele: reads with its base, coverage, mean base quality."""
    names = list(aln)
    L = len(aln[names[0]])
    pos = {n: [] for n in names}                            # column -> 0-based allele coordinate
    for n in names:
        i = -1
        for c in aln[n]:
            i += c != "-"
            pos[n].append(i if c != "-" else None)
    var = [c for c in range(L) if len({aln[n][c] for n in names}) > 1]
    rows = []
    for n in names:
        bam = pysam.AlignmentFile(bams[n])
        want = {pos[n][c]: c for c in var if pos[n][c] is not None}
        for col in bam.pileup(n, truncate=True, min_base_quality=0, stepper="nofilter"):
            if col.reference_pos not in want:
                continue
            c = want[col.reference_pos]
            base = aln[n][c]
            cov = hit = 0; q = []
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip:
                    continue
                if keep is not None and pr.alignment.query_name not in keep:
                    continue   # only reads whose best hit is this gene (no paralog cross-mapping)
                cov += 1
                b = pr.alignment.query_sequence[pr.query_position]
                if b == base:
                    hit += 1; q.append(pr.alignment.query_qualities[pr.query_position])
            rows.append({"allele": n, "column": c + 1, "allele_pos": col.reference_pos + 1, "base": base,
                         "coverage": cov, "reads_with_base": hit,
                         "mean_q_with_base": round(sum(q) / len(q), 1) if q else None})
        for c in var:                                       # sites without coverage
            if pos[n][c] is None or not any(r["allele"] == n and r["column"] == c + 1 for r in rows):
                rows.append({"allele": n, "column": c + 1, "allele_pos": (pos[n][c] or -1) + 1,
                             "base": aln[n][c], "coverage": 0, "reads_with_base": 0, "mean_q_with_base": None})
    return pd.DataFrame(rows), var


def plot(sample, genes, aln_by_gene, sites, sup, calls, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not genes:
        return
    fig, axes = plt.subplots(len(genes), 2, figsize=(13, 3.4 * len(genes)), squeeze=False,
                             gridspec_kw={"width_ratios": [3, 1]})
    for ax, g in zip(axes, genes):
        s = sites[sites.gene == g]
        cols = sorted(s.column.unique())
        names = list(aln_by_gene[g])
        for i, n in enumerate(names):
            d = s[s.allele == n].set_index("column").reindex(cols)
            frac = (d.reads_with_base / d.coverage.where(d.coverage > 0))
            ax[0].plot(range(len(cols)), frac.values, "o-", ms=3, lw=.8, label=n)
        ax[0].set_xticks(range(len(cols)))
        ax[0].set_xticklabels(cols, rotation=90, fontsize=6)
        ax[0].set_ylim(-.05, 1.05)
        ax[0].set_ylabel("reads with allele base / coverage", fontsize=8)
        ax[0].set_xlabel("variable MAFFT column (CDS)", fontsize=8)
        ax[0].legend(fontsize=7, frameon=False, ncol=2)
        ax[0].set_title(f"{sample} HLA-{g}: per-site read support", fontsize=9)
        t = sup[sup.gene == g].sort_values("explained", ascending=False).head(6)
        lab = [f"{r.allele1}/{r.allele2}".replace(f"{g}*", "") for _, r in t.iterrows()]
        col = ["#1f77b4" if [r.allele1, r.allele2] == calls[g] else "#bbbbbb" for _, r in t.iterrows()]
        ax[1].barh(lab[::-1], t.pct_explained.values[::-1], color=col[::-1])
        ax[1].set_xlabel("% gene reads explained", fontsize=8)
        ax[1].tick_params(labelsize=7)
        ax[1].set_title("genotype support (blue = final)", fontsize=9)
    plt.tight_layout(); plt.savefig(out, dpi=150); plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--optitype", required=True, help="OptiType *_result.tsv")
    ap.add_argument("--arcas", required=True, help="arcasHLA *.genotype.json")
    ap.add_argument("--reads", nargs="+", required=True, help="HLA reads (fished FASTQ, R1 [R2])")
    ap.add_argument("--hla-nuc", required=True)
    ap.add_argument("--cellosaurus", default=None, help="CVCL_ accession or exact cell-line name")
    ap.add_argument("--cellosaurus-file", default=None, help="user-typed TSV: gene allele1 allele2 [source]")
    ap.add_argument("--override", default=None, help="user-typed final calls TSV: gene allele1 allele2")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    od = Path(a.outdir); od.mkdir(parents=True, exist_ok=True)
    s = a.sample

    opt, opt_reads = optitype_calls(a.optitype)
    arc = arcas_calls(a.arcas)
    cac, crecs = cellosaurus_records(a)
    nuc = read_nuc(a.hla_nuc, set(GENES))

    cand_by_gene, cands, skipped = defaultdict(set), {}, []
    for g in GENES:
        pool = set(opt[g]) | set(arc.get(g, []))
        for rec, _ in crecs:
            pool |= set(rec.get(g, []))
        for t in pool:
            if t.count(":") < 1:
                continue                                    # 1-field Cellosaurus entries: comparison only
            try:
                full = pick(t, nuc)
            except SystemExit as e:
                skipped.append(str(e)); continue
            cands[t] = nuc[full]; cand_by_gene[g].add(t)

    with tempfile.TemporaryDirectory() as tmp:
        reads = f"{tmp}/reads.fq"
        with open(reads, "w") as fo:                        # mates as single reads, unique names
            for i, r in enumerate(a.reads):
                for rec in pysam.FastxFile(r):
                    fo.write(f"@{rec.name}/{i + 1}\n{rec.sequence}\n+\n{rec.quality}\n")
        bams = align_each(cands, reads, tmp, a.threads)
        sup, gene_reads = genotype_support(best_as(bams), cand_by_gene)

        final, rows, dis_genes, aln_by_gene, site_tabs = {}, [], [], {}, []
        over = pd.read_csv(a.override, sep="\t").set_index("gene") if a.override else None
        for g in GENES:
            sg = sup[sup.gene == g].set_index(["allele1", "allele2"])
            ex = lambda p: int(sg.loc[tuple(p), "explained"]) if p and tuple(p) in sg.index else None
            o, r_ = opt[g], arc.get(g)
            cel = [(rec[g], src) for rec, src in crecs if g in rec]
            cel2 = sorted([(tuple(c), src, ex(c)) for c, src in cel if all(":" in x for x in c)],
                          key=lambda t: -(t[2] or -1))       # best read-supported record first
            cbest = list(cel2[0][0]) if cel2 else None

            # site evidence: alleles of the disputed pairs, at columns where they differ
            disputed = sorted(set(o) ^ set(r_ or []) | (set(cbest) ^ set(o) if cbest else set()))
            disputed = [t for t in disputed if t in cands]
            rej = {}
            if disputed:
                aln = mafft({t: cands[t][1] for t in sorted(cand_by_gene[g])}, tmp, g)
                st, _ = site_evidence(aln, bams, keep=set(gene_reads[g]))
                cols = {c for c in st.column.unique()
                        if len({aln[t][c - 1] for t in disputed + list(set(o) & set(r_ or []))}) > 1}
                st = st[st.column.isin(cols)].copy()
                st["frac"] = st.reads_with_base / st.coverage.where(st.coverage > 0)
                st.insert(0, "gene", g); site_tabs.append(st); aln_by_gene[g] = aln; dis_genes.append(g)
                for t in disputed:                           # allele-specific base barely seen -> reject
                    d = st[st.allele == t]
                    shared = d.column.map(lambda c: sum(aln[u][c - 1] == aln[t][c - 1] for u in aln) > 1)
                    d = d[~shared.values] if (~shared).any() else d
                    f = d.frac.median()
                    worst = d[d.coverage >= MIN_SITE_COV].frac.min()
                    if (pd.notna(f) and f < MIN_SITE_FRAC) or (pd.notna(worst) and worst < MIN_SITE_FRAC_ANY):
                        rej[t] = f"median {f:.3f}, min {worst:.3f}"

            ok = lambda p: p is not None and not (set(p) & set(rej))
            tie = max(MIN_DELTA_READS, MIN_DELTA_PCT / 100 * int(sg.gene_reads.iloc[0]))
            if over is not None and g in over.index:
                call, why = sorted([two(over.loc[g, "allele1"]), two(over.loc[g, "allele2"])]), "user_override"
            elif o == r_ and ok(o):
                call, why = o, "tools_agree"
            else:
                pool = [(p, n) for p, n in [(o, "OptiType"), (r_, "arcasHLA"), (cbest, "Cellosaurus")] if ok(p)]
                if not pool:
                    call, why = o, "REVIEW_all_rejected"
                else:
                    top = max(ex(p) or 0 for p, _ in pool)
                    near = [(p, n) for p, n in pool if top - (ex(p) or 0) < tie]
                    if len({tuple(p) for p, _ in near}) == 1:
                        call, why = near[0][0], f"{near[0][1]}_higher_support"
                    else:                                    # near tie: Cellosaurus > OptiType > arcasHLA
                        call, n = next(((p, n) for p, n in near if n == "Cellosaurus"), near[0])
                        why = f"near_tie_{n}_REVIEW"
            final[g] = call
            best = sup[sup.gene == g].sort_values(["explained", "AS_lost"], ascending=[False, True]).iloc[0]
            rows.append({"sample": s, "gene": g, "allele1": call[0], "allele2": call[1], "decision": why,
                         "optitype": "/".join(o), "arcasHLA": "/".join(r_ or []),
                         "support_optitype": ex(o), "support_arcasHLA": ex(r_), "support_final": ex(call),
                         "gene_reads": int(best.gene_reads),
                         "best_supported": f"{best.allele1}/{best.allele2}", "support_best": int(best.explained),
                         "rejected_by_sites": "; ".join(f"{k} ({v})" for k, v in rej.items()),
                         "cellosaurus_best": "/".join(cbest) if cbest else "",
                         "cellosaurus_best_source": cel2[0][1] if cel2 else "",
                         "cellosaurus_match": (call == cbest) if cbest else None,
                         "cellosaurus_records": "; ".join(f"{'/'.join(c)} [{src}] n={n}" for c, src, n in cel2)})

    calls = pd.DataFrame(rows)
    calls.to_csv(od / f"{s}_hla_calls.tsv", sep="\t", index=False)
    sup.insert(0, "sample", s); sup.to_csv(od / f"{s}_hla_support.tsv", sep="\t", index=False)
    sites = pd.concat(site_tabs) if site_tabs else pd.DataFrame(columns=["gene", "allele", "column"])
    sites.to_csv(od / f"{s}_hla_sites.tsv", sep="\t", index=False)
    pd.DataFrame([{"accession": cac, "source": src, **{g: "/".join(rec.get(g, [])) for g in GENES}}
                  for rec, src in crecs]).to_csv(od / f"{s}_hla_cellosaurus.tsv", sep="\t", index=False)
    plot(s, dis_genes, aln_by_gene, sites, sup, final, od / f"{s}_hla_disagreement.png")
    for m in skipped:
        print("skipped candidate:", m, file=sys.stderr)
    print(calls.drop(columns=["sample"]).to_string(index=False))


if __name__ == "__main__":
    main()
