#!/usr/bin/env python3
"""Classify genes by WHY they are ambiguous and WHAT to do about it.

Two independent questions, answered separately:

  cause   - from annotation geometry and identical terminal sequence:
              alt_haplotype_copy    shares it with another copy of the SAME gene name, i.e. an
                                    alt haplotype, patch or unplaced scaffold in the annotation
              pseudogene_shadow     shares terminal sequence with a pseudogene
              local_tandem_dup      shares it with a coding gene, same chromosome, within --local-mb
              distant_paralog       shares it with a coding gene, far away or another chromosome
              readthrough           annotation tag
              retrogene             annotation tag
              nested                interval inside a longer gene
              overlap_antisense     shares coordinates, opposite strand
              overlap_sense         shares coordinates, same strand
              self_homology         internal repeat (set via --self-homology)
              unique                none of the above
            A gene can have several; `cause` is the highest-priority one and `all_causes` lists all.

  action  - from funique and, when given, the inflation ratio between a UTR-containing and a
            CDS-only reference:
              always_ok            funique >= --uf-high, no inflation      any pipeline is fine
              apply_correction     funique below the cut, no inflation     multiply the count by funique
              prefer_pseudoaligner funique high but inflated               unique counting is right
              always_wrong         funique <= --uf-zero and inflated       neither works
            Correction is advised ONLY when the competing sequence is outside the reference, i.e.
            for a pseudogene cause. When the competitor is inside the reference its reads are
            already divided by fractional counting, and scaling again corrects twice.

usage: classify_genes.py --geometry F --clusters F --funique F --tx2gene F --out F
                         [--inflation F] [--self-homology F] [--local-mb 5]
"""
import argparse

import pandas as pd

PRIORITY = ["alt_haplotype_copy", "pseudogene_shadow", "local_tandem_dup", "distant_paralog",
            "self_homology", "retrogene", "readthrough", "nested", "overlap_sense",
            "overlap_antisense", "unique"]
# causes where the competing sequence is absent from a one-transcript-per-gene reference
EXTERNAL = {"pseudogene_shadow", "retrogene"}
# causes that mean the reference itself is wrong: fix the reference, do not scale counts
REFERENCE_FAULT = {"alt_haplotype_copy"}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--geometry", required=True)
    p.add_argument("--clusters", required=True)
    p.add_argument("--funique", required=True)
    p.add_argument("--tx2gene", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--inflation", help="TSV: gene_id, inflation_ratio (UTR ref / CDS ref)")
    p.add_argument("--self-homology", help="TSV with gene_id of genes carrying internal repeats")
    p.add_argument("--coding-type", default="protein_coding")
    p.add_argument("--local-mb", type=float, default=5.0)
    p.add_argument("--uf-high", type=float, default=0.9)
    p.add_argument("--uf-zero", type=float, default=0.1)
    p.add_argument("--inflation-hi", type=float, default=2.0)
    a = p.parse_args()

    geo = pd.read_csv(a.geometry, sep="\t")
    cl = pd.read_csv(a.clusters, sep="\t")
    uf = pd.read_csv(a.funique, sep="\t")
    t2g = pd.read_csv(a.tx2gene, sep="\t")

    # funique is per transcript; carry it to the gene via the representative transcript
    uf = uf.merge(t2g[["transcript_id", "gene_id"]], on="transcript_id", how="left")
    uf_gene = uf.groupby("gene_id", as_index=False).agg(
        funique=("funique", "max"), n_positions=("n_positions", "max"))

    pos = geo.set_index("gene_id")[["chrom", "start", "end", "gene_type"]].to_dict("index")
    is_pseudo = lambda t: "pseudogene" in str(t)

    # cluster partners -> cause
    cause_sets = {g: set() for g in geo.gene_id}
    partner_note = {}
    for cid, grp in cl.groupby("cluster_id"):
        gids = grp.gene_id.unique()
        if len(gids) < 2:
            continue
        for g in gids:
            if g not in cause_sets:
                continue
            others = [o for o in gids if o != g]
            gi = pos.get(g)
            gname = grp.loc[grp.gene_id == g, "gene_name"].iloc[0]
            for o in others:
                oi = pos.get(o)
                ot = (oi or {}).get("gene_type", "") or grp.loc[grp.gene_id == o, "biotype"].iloc[0]
                oname = grp.loc[grp.gene_id == o, "gene_name"].iloc[0]
                if oname == gname:
                    # same symbol, different gene_id: a second copy of this gene in the
                    # annotation (alt haplotype / patch / scaffold), not a paralog
                    cause_sets[g].add("alt_haplotype_copy")
                elif is_pseudo(ot):
                    cause_sets[g].add("pseudogene_shadow")
                elif gi and oi and oi["chrom"] == gi["chrom"] and \
                        abs(oi["start"] - gi["start"]) <= a.local_mb * 1e6:
                    cause_sets[g].add("local_tandem_dup")
                else:
                    cause_sets[g].add("distant_paralog")
            partner_note[g] = ",".join(sorted(
                set(grp.loc[grp.gene_id.isin(others), "gene_name"].astype(str)))[:20])

    if a.self_homology:
        sh = pd.read_csv(a.self_homology, sep="\t")
        for g in sh.gene_id:
            if g in cause_sets:
                cause_sets[g].add("self_homology")

    for r in geo.itertuples():
        tg = str(r.tags)
        if "readthrough" in tg:
            cause_sets[r.gene_id].add("readthrough")
        if "retrogene" in tg:
            cause_sets[r.gene_id].add("retrogene")
        if r.is_nested:
            cause_sets[r.gene_id].add("nested")
        if r.overlap_group_size > 1:
            cause_sets[r.gene_id].add("overlap_antisense" if r.overlap_strand_config == "both_strands"
                                      else "overlap_sense")
        if not cause_sets[r.gene_id]:
            cause_sets[r.gene_id].add("unique")

    top = {g: next(c for c in PRIORITY if c in s) for g, s in cause_sets.items()}

    d = geo.merge(uf_gene, on="gene_id", how="left")
    d["cause"] = d.gene_id.map(top)
    d["all_causes"] = d.gene_id.map(lambda g: ";".join(sorted(cause_sets.get(g, {"unique"}))))
    d["cluster_partners"] = d.gene_id.map(partner_note).fillna("")

    if a.inflation:
        inf = pd.read_csv(a.inflation, sep="\t")
        d = d.merge(inf[["gene_id", "inflation_ratio"]], on="gene_id", how="left")
    else:
        d["inflation_ratio"] = pd.NA

    def action(r):
        f = r.funique
        if pd.isna(f):
            return "no_data"
        inflated = pd.notna(r.inflation_ratio) and r.inflation_ratio >= a.inflation_hi
        external = r.cause in EXTERNAL
        if r.cause in REFERENCE_FAULT:
            return "use_primary_assembly"
        if f <= a.uf_zero:
            return "always_wrong" if inflated else "unquantifiable"
        if f >= a.uf_high:
            return "prefer_pseudoaligner" if inflated else "always_ok"
        return "apply_correction" if external else "ambiguous_no_correction"

    d["action"] = d.apply(action, axis=1)
    d.loc[d.gene_type != a.coding_type, "action"] = "not_coding"

    cols = ["gene_id", "gene_name", "gene_type", "chrom", "start", "end", "strand",
            "cause", "all_causes", "cluster_partners", "overlap_group_size",
            "overlap_n_pseudogene", "is_nested", "funique", "inflation_ratio", "action"]
    d[cols].to_csv(a.out, sep="\t", index=False)

    coding = d[d.gene_type == a.coding_type]
    print(f"{len(d)} genes, {len(coding)} {a.coding_type}\n")
    print("cause (coding genes):")
    print(coding.cause.value_counts().to_string())
    print("\nmean funique by cause:")
    print(coding.groupby("cause").funique.agg(["size", "mean"]).round(3).to_string())
    print("\naction:")
    print(coding.action.value_counts().to_string())


if __name__ == "__main__":
    main()
