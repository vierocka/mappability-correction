#!/usr/bin/env python3
"""
Compare, per gene symbol, how many distinct CDS-bearing transcripts each of
three annotation sources assigns to that gene:

  - Homo_sapiens.GRCh38.116.gtf        (Ensembl, full annotation, all isoforms)
  - gencode.v49.basic.annotation.gtf   (GENCODE "basic" tag — curated subset)
  - GCF_000001405.40_GRCh38.p14_genomic.gff  (NCBI RefSeq, GRCh38.p14)

Counts transcripts with >=1 CDS feature, not individual CDS exon lines — i.e.
how many coding-transcript models each source annotates for the gene, which
is what actually drives dedup-CDS reference size (see extract_dedup.log).

Ensembl/GENCODE use GTF format (key "value"; attributes); RefSeq uses GFF3
(key=value attributes, CDS Parent= points at an rna- feature ID rather than
repeating transcript_id directly) — parsed with two different attribute
parsers accordingly.

Output: one TSV, one row per gene symbol seen in any of the three sources.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_gtf_attrs(col9):
    """Ensembl/GENCODE GTF: key "value"; ... (space-separated key/value pairs)."""
    attrs = {}
    for field in col9.strip().split(';'):
        field = field.strip()
        if not field:
            continue
        key, _, val = field.partition(' ')
        attrs[key] = val.strip('"')
    return attrs


def parse_gff3_attrs(col9):
    """NCBI GFF3: key=value;key=value (no quotes)."""
    attrs = {}
    for field in col9.strip().split(';'):
        if '=' in field:
            k, v = field.split('=', 1)
            attrs[k] = v
    return attrs


def count_gtf(path):
    """{gene_symbol: n_distinct_CDS-bearing_transcripts} for an Ensembl/GENCODE-style GTF."""
    gene_tx = defaultdict(set)
    n_lines = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'CDS':
                continue
            a = parse_gtf_attrs(f[8])
            gene = a.get('gene_name')
            tid = a.get('transcript_id')
            if gene is None or tid is None:
                continue
            gene_tx[gene].add(tid)
            n_lines += 1
    print(f"    CDS lines: {n_lines:,}  genes with CDS: {len(gene_tx):,}", flush=True)
    return {g: len(txs) for g, txs in gene_tx.items()}


def count_gff3(path):
    """{gene_symbol: n_distinct_CDS-bearing_transcripts} for the NCBI RefSeq GFF3."""
    gene_tx = defaultdict(set)
    n_lines = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'CDS':
                continue
            a = parse_gff3_attrs(f[8])
            gene = a.get('gene')
            tid = a.get('Parent')  # e.g. rna-NM_000123.4 — one CDS can span multiple lines (exons), same Parent
            if gene is None or tid is None:
                continue
            gene_tx[gene].add(tid)
            n_lines += 1
    print(f"    CDS lines: {n_lines:,}  genes with CDS: {len(gene_tx):,}", flush=True)
    return {g: len(txs) for g, txs in gene_tx.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ensembl-gtf', help='Homo_sapiens.GRCh38.116.gtf')
    ap.add_argument('--gencode-gtf', help='gencode.v49.basic.annotation.gtf')
    ap.add_argument('--refseq-gff', help='GCF_000001405.40_GRCh38.p14_genomic.gff')
    ap.add_argument('--out', default='cds_counts_by_gene.tsv')
    args = ap.parse_args()

    if not any([args.ensembl_gtf, args.gencode_gtf, args.refseq_gff]):
        ap.error('provide at least one of --ensembl-gtf / --gencode-gtf / --refseq-gff')

    ens, genc, refseq = {}, {}, {}

    if args.ensembl_gtf:
        print(f"Parsing Ensembl GTF: {args.ensembl_gtf}", flush=True)
        ens = count_gtf(Path(args.ensembl_gtf))
    else:
        print("Skipping Ensembl GTF (not provided)")

    if args.gencode_gtf:
        print(f"Parsing GENCODE basic GTF: {args.gencode_gtf}", flush=True)
        genc = count_gtf(Path(args.gencode_gtf))
    else:
        print("Skipping GENCODE GTF (not provided)")

    if args.refseq_gff:
        print(f"Parsing RefSeq p14 GFF: {args.refseq_gff}", flush=True)
        refseq = count_gff3(Path(args.refseq_gff))
    else:
        print("Skipping RefSeq GFF (not provided)")

    all_genes = sorted(set(ens) | set(genc) | set(refseq))
    with open(args.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['gene_symbol', 'n_cds_transcripts_ensembl116',
                    'n_cds_transcripts_gencode_v49_basic', 'n_cds_transcripts_refseq_p14'])
        for g in all_genes:
            w.writerow([g, ens.get(g, 0), genc.get(g, 0), refseq.get(g, 0)])

    print(f"\nWrote {len(all_genes):,} genes -> {args.out}")


if __name__ == '__main__':
    main()
