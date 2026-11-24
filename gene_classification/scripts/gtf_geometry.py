#!/usr/bin/env python3
"""Gene geometry from an annotation GTF: overlap groups, nesting, biotype, tags.

Organism-agnostic. Reads any GTF that carries gene records with a gene id and a biotype
attribute (`gene_type` as in GENCODE, or `gene_biotype` as in Ensembl/RefSeq-derived GTFs).

Produces one row per gene:
  gene_id gene_name gene_type chrom start end strand tags
  overlap_group_id overlap_group_size overlap_strand_config
  overlap_n_coding overlap_n_pseudogene overlap_group_genes is_nested

An overlap group is a connected component of genes sharing genomic coordinates, found by an
interval sweep per chromosome. is_nested is set when a gene's interval lies entirely inside
another gene's interval.

usage: gtf_geometry.py --gtf FILE --out FILE [--coding-type protein_coding]
"""
import argparse
import gzip
import re
from collections import defaultdict

OPEN = lambda p: gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def attr(s, key):
    m = re.search(key + r'\s+"([^"]*)"', s)
    return m.group(1) if m else ""


def all_tags(s):
    return ";".join(re.findall(r'tag\s+"([^"]*)"', s))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gtf", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--coding-type", default="protein_coding")
    a = p.parse_args()

    genes = []
    with OPEN(a.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            gt = attr(c[8], "gene_type") or attr(c[8], "gene_biotype")
            gid = attr(c[8], "gene_id")
            if not gid:
                continue
            genes.append({
                "gene_id": gid,
                "gene_name": attr(c[8], "gene_name") or gid,
                "gene_type": gt,
                "chrom": c[0], "start": int(c[3]), "end": int(c[4]), "strand": c[6],
                "tags": all_tags(c[8]),
            })

    if not genes:
        raise SystemExit("no gene records found - does this GTF have 'gene' features?")

    by_chrom = defaultdict(list)
    for g in genes:
        by_chrom[g["chrom"]].append(g)

    # interval sweep -> connected components of overlapping intervals
    gid_group, group_members = {}, defaultdict(list)
    gcount = 0
    for chrom, gl in by_chrom.items():
        gl.sort(key=lambda g: (g["start"], g["end"]))
        cur, cur_end = [], -1
        for g in gl:
            if cur and g["start"] > cur_end:          # gap: close the component
                gcount += 1
                for m in cur:
                    gid_group[m["gene_id"]] = f"G{gcount}"
                    group_members[f"G{gcount}"].append(m)
                cur, cur_end = [], -1
            cur.append(g)
            cur_end = max(cur_end, g["end"])
        if cur:
            gcount += 1
            for m in cur:
                gid_group[m["gene_id"]] = f"G{gcount}"
                group_members[f"G{gcount}"].append(m)

    is_pseudo = lambda t: "pseudogene" in t

    rows = []
    for g in genes:
        grp = gid_group[g["gene_id"]]
        members = group_members[grp]
        others = [m for m in members if m["gene_id"] != g["gene_id"]]
        strands = {m["strand"] for m in members}
        nested = any(m["start"] <= g["start"] and m["end"] >= g["end"]
                     and (m["end"] - m["start"]) > (g["end"] - g["start"]) for m in others)
        rows.append({
            **g,
            "overlap_group_id": grp,
            "overlap_group_size": len(members),
            "overlap_strand_config": ("single" if len(members) == 1 else
                                      "same_strand" if len(strands) == 1 else "both_strands"),
            "overlap_n_coding": sum(m["gene_type"] == a.coding_type for m in members),
            "overlap_n_pseudogene": sum(is_pseudo(m["gene_type"]) for m in members),
            "overlap_group_genes": ",".join(sorted(m["gene_name"] for m in others)[:20]),
            "is_nested": int(nested),
        })

    cols = ["gene_id", "gene_name", "gene_type", "chrom", "start", "end", "strand", "tags",
            "overlap_group_id", "overlap_group_size", "overlap_strand_config",
            "overlap_n_coding", "overlap_n_pseudogene", "overlap_group_genes", "is_nested"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    n_multi = sum(r["overlap_group_size"] > 1 for r in rows)
    print(f"{len(rows)} genes, {gcount} overlap groups, "
          f"{n_multi} genes ({100*n_multi/len(rows):.1f} %) share coordinates with another gene, "
          f"{sum(r['is_nested'] for r in rows)} nested")


if __name__ == "__main__":
    main()
