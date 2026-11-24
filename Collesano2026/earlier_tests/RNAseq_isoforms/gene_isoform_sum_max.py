#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass

@dataclass
class IsoRec:
    isoform: str
    length: int
    est_counts: float

def parse_args():
    ap = argparse.ArgumentParser(
        description="Aggregate kallisto isoform est_counts to gene-level sum and max (tie -> longest)."
    )
    ap.add_argument("--abundance", required=True, help="kallisto abundance.tsv")
    ap.add_argument("--annot", required=True, help="Tab-separated annotation: CDS<tab>gene<tab>xref")
    ap.add_argument("--out", required=True, help="Output TSV")
    ap.add_argument("--use", choices=["est_counts", "tpm"], default="est_counts",
                    help="Which abundance to aggregate (default: est_counts)")
    return ap.parse_args()

def extract_isoform_id(target_id: str) -> str | None:
    """
    target_id example:
      lcl|NC_000001.11_cds_NP_001005484.2_1
    Extract isoform ID between '_cds_' and the last underscore-suffix.
    """
    if "_cds_" not in target_id:
        return None
    after = target_id.split("_cds_", 1)[1]          # e.g. NP_001005484.2_1
    parts = after.split("_")
    if len(parts) < 2:
        return None
    isoform = "_".join(parts[:-1])                 # drop last numeric suffix
    return isoform

def read_annotation_map(annot_path: str) -> dict[str, str]:
    """
    Reads tab-separated file with header:
      CDS    gene    xref
    Returns isoform(CDS) -> gene_name
    """
    iso2gene: dict[str, str] = {}
    with open(annot_path, "rt", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None or "CDS" not in reader.fieldnames or "gene" not in reader.fieldnames:
            raise ValueError("Annotation file must be tab-separated with columns: CDS, gene, xref")
        for row in reader:
            cds = (row.get("CDS") or "").strip()
            gene = (row.get("gene") or "").strip()
            if cds and gene and cds not in iso2gene:
                iso2gene[cds] = gene  # first occurrence wins
    return iso2gene

def main():
    args = parse_args()
    iso2gene = read_annotation_map(args.annot)

    # gene -> list of isoform records (we’ll fold on the fly)
    gene_sum = defaultdict(float)

    # gene -> current best isoform record according to:
    # (est_counts, length) with est_counts primary, length tie-break
    gene_best: dict[str, IsoRec] = {}

    use_col = args.use  # "est_counts" or "tpm"

    with open(args.abundance, "rt", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        needed = {"target_id", "length", use_col}
        if reader.fieldnames is None or not needed.issubset(set(reader.fieldnames)):
            raise ValueError(f"abundance.tsv must contain columns: {sorted(needed)}")

        for row in reader:
            target_id = (row.get("target_id") or "").strip()
            iso = extract_isoform_id(target_id)
            if not iso:
                continue

            gene = iso2gene.get(iso)
            if gene is None:
                continue  # skip isoforms without gene mapping

            # parse numbers
            try:
                length = int(float(row["length"]))
            except Exception:
                length = 0

            try:
                val = float(row[use_col])
            except Exception:
                continue

            gene_sum[gene] += val

            rec = IsoRec(isoform=iso, length=length, est_counts=val)

            if gene not in gene_best:
                gene_best[gene] = rec
            else:
                best = gene_best[gene]
                # choose higher val; tie -> longer length; tie -> keep first (do nothing)
                if rec.est_counts > best.est_counts:
                    gene_best[gene] = rec
                elif rec.est_counts == best.est_counts and rec.length > best.length:
                    gene_best[gene] = rec

    # write output
    with open(args.out, "wt", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow([
            "gene",
            f"sum_{use_col}",
            f"max_{use_col}",
            "best_isoform",
            "best_isoform_length",
        ])
        for gene in sorted(gene_sum.keys()):
            best = gene_best.get(gene)
            if best is None:
                continue
            w.writerow([
                gene,
                f"{gene_sum[gene]:.6f}",
                f"{best.est_counts:.6f}",
                best.isoform,
                str(best.length),
            ])

if __name__ == "__main__":
    main()
