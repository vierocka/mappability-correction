#!/usr/bin/env python3
"""9-mer table (gene_id, gene_symbol, peptide_9mer, count_in_gene) -> one CSV row per
occurrence (transcript_id, gene_symbol, peptide): the row layout rebuild_peptide_tables.py sums."""
import csv, sys
src, dst = sys.argv[1], sys.argv[2]
with open(src) as f, open(dst, "w", newline="") as o:
    r, w = csv.DictReader(f, delimiter="\t"), csv.writer(o)
    w.writerow(["transcript_id", "gene_symbol", "peptide"])
    for x in r:
        for _ in range(int(x["count_in_gene"])):
            w.writerow([x["gene_id"], x["gene_symbol"], x["peptide_9mer"]])
