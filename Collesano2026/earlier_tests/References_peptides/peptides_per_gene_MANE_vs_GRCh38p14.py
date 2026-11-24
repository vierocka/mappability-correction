#!/usr/bin/env python3
"""
peptides_per_gene.py
====================
Compare unique 9-mer peptides per gene between:
  1. MANE.GRCh38.v1.5.ensembl_protein.faa        (one canonical isoform per gene)
  2. GCF_000001405.40_GRCh38.p14_protein.dedup.faa (multiple isoforms per gene)

The GCF set requires CDS_GRCh38p14_annotNames.csv (tab-separated, columns: CDS, gene, xref)
to map protein IDs to gene symbols.

Outputs (TSV):
  mane_9mers.tsv   — per-gene summary for the MANE set
  gcf_9mers.tsv    — per-gene summary for the GCF set

Both output files share the same columns:
  geneID | num_isoforms | num_unique_9mers | unique_9mers (comma-separated)

Note: the unique_9mers column can be very large. If you only need counts,
set WRITE_PEPTIDE_LIST = False below.

Usage:
  python peptides_per_gene.py
  (all input files must be in the current working directory)

Author:  generated for the immunopeptidomics pipeline
"""

import sys
import csv
from collections import defaultdict

# ── Configuration ────────────────────────────────────────────────────────────

MANE_FAA   = "<SELF_NONSELF>/human/References/MANE_v1.5/MANE.GRCh38.v1.5.ensembl_protein.faa"
GCF_FAA    = "<SELF_NONSELF>/human/References/GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.dedup.faa"
ANNOT_CSV  = "<SELF_NONSELF>/human/References/GRCh38.p14/CDS_GRCh38p14_annotNames.csv"   # tab-separated

KMER_LEN   = 9          # peptide length (change here if you ever need 8- or 10-mers)
WRITE_PEPTIDE_LIST = True   # set False to skip the (potentially huge) peptide column

# ── Helpers ───────────────────────────────────────────────────────────────────

def parse_fasta(filepath):
    """Yield (header_string, sequence) for every record in a FASTA file."""
    header, parts = None, []
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts)
                header, parts = line[1:], []
            else:
                parts.append(line)
    if header is not None:
        yield header, ''.join(parts)


def get_kmers(sequence, k=KMER_LEN):
    """Return a set of all k-mer substrings from sequence."""
    n = len(sequence)
    if n < k:
        return set()
    return {sequence[i:i+k] for i in range(n - k + 1)}


def parse_mane_gene(header):
    """
    Extract gene symbol from a MANE / Ensembl protein FASTA header.

    Typical format:
      ENSP00000493376.2 pep chromosome:GRCh38:... gene:ENSG00000186092.6 ... gene_symbol:OR4F5 ...

    Priority: gene_symbol: > gene: > first token (protein ID as last resort).
    """
    gene_symbol = None
    gene_id     = None
    for token in header.split():
        if token.startswith('gene_symbol:'):
            gene_symbol = token.split(':', 1)[1]
        elif token.startswith('gene:'):
            gene_id = token.split(':', 1)[1]
    return gene_symbol or gene_id or header.split()[0]


def load_cds_gene_map(csv_path):
    """
    Parse CDS_GRCh38p14_annotNames.csv (tab-separated, header row).
    Returns dict: protein_id (e.g. NP_001005484.2) -> gene symbol (e.g. OR4F5).
    """
    mapping = {}
    with open(csv_path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            cds  = row['CDS'].strip()
            gene = row['gene'].strip()
            if cds and gene:
                mapping[cds] = gene
    return mapping


def process_fasta(fasta_path, gene_fn):
    """
    Read a FASTA file, group sequences by gene (using gene_fn to extract gene from header),
    and accumulate unique k-mers per gene.

    gene_fn: callable(header) -> gene_id string  (or None to skip the record)

    Returns dict: gene -> {'isoforms': int, 'peptides': set}
    """
    genes    = defaultdict(lambda: {'isoforms': 0, 'peptides': set()})
    skipped  = 0
    total    = 0

    for header, seq in parse_fasta(fasta_path):
        total += 1
        gene = gene_fn(header)
        if gene is None:
            skipped += 1
            continue
        genes[gene]['isoforms'] += 1
        genes[gene]['peptides'].update(get_kmers(seq))

    return genes, total, skipped


def write_output(outpath, gene_data):
    """Write TSV output for one protein set."""
    cols = ['geneID', 'num_isoforms', 'num_unique_9mers']
    if WRITE_PEPTIDE_LIST:
        cols.append('unique_9mers')

    with open(outpath, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for gene in sorted(gene_data):
            d = gene_data[gene]
            peptides = sorted(d['peptides'])
            row = [gene, str(d['isoforms']), str(len(peptides))]
            if WRITE_PEPTIDE_LIST:
                row.append(','.join(peptides))
            fh.write('\t'.join(row) + '\n')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():

    # ── 1. MANE set ──────────────────────────────────────────────────────────
    print(f"[1/2] Processing MANE set: {MANE_FAA}", file=sys.stderr)

    mane_data, mane_total, mane_skip = process_fasta(MANE_FAA, parse_mane_gene)

    write_output("mane_9mers.tsv", mane_data)

    print(f"      {mane_total} proteins read | "
          f"{mane_skip} skipped | "
          f"{len(mane_data)} genes → mane_9mers.tsv", file=sys.stderr)

    # Quick sanity: MANE should be 1 isoform per gene
    multi = [(g, d['isoforms']) for g, d in mane_data.items() if d['isoforms'] > 1]
    if multi:
        print(f"      WARNING: {len(multi)} genes have >1 isoform in MANE "
              f"(e.g. {multi[:3]})", file=sys.stderr)

    # ── 2. GCF deduplicated set ───────────────────────────────────────────────
    print(f"[2/2] Processing GCF set:  {GCF_FAA}", file=sys.stderr)
    print(f"      Loading gene map from {ANNOT_CSV} …", file=sys.stderr)

    cds_to_gene = load_cds_gene_map(ANNOT_CSV)
    print(f"      {len(cds_to_gene)} CDS → gene mappings loaded", file=sys.stderr)

    def gcf_gene_fn(header):
        """Extract protein ID from GCF header and look up gene in the mapping."""
        protein_id = header.split()[0]   # e.g.  NP_001005484.2
        return cds_to_gene.get(protein_id)  # None if not found → record skipped

    gcf_data, gcf_total, gcf_skip = process_fasta(GCF_FAA, gcf_gene_fn)

    write_output("gcf_9mers.tsv", gcf_data)

    print(f"      {gcf_total} proteins read | "
          f"{gcf_skip} unmapped/skipped | "
          f"{len(gcf_data)} genes → gcf_9mers.tsv", file=sys.stderr)
    if gcf_skip:
        print(f"      TIP: unmapped proteins are absent from {ANNOT_CSV}. "
              f"This may include proteins added after the annotation snapshot.",
              file=sys.stderr)

    # ── 3. Brief cross-set summary ────────────────────────────────────────────
    print("\n── Summary ──────────────────────────────────────────────", file=sys.stderr)
    shared = set(mane_data) & set(gcf_data)
    mane_only = set(mane_data) - set(gcf_data)
    gcf_only  = set(gcf_data)  - set(mane_data)

    mane_kmer_total = sum(len(d['peptides']) for d in mane_data.values())
    gcf_kmer_total  = sum(len(d['peptides']) for d in gcf_data.values())

    print(f"  Genes in MANE only  : {len(mane_only):>6}", file=sys.stderr)
    print(f"  Genes in GCF only   : {len(gcf_only):>6}", file=sys.stderr)
    print(f"  Genes in both sets  : {len(shared):>6}", file=sys.stderr)
    print(f"  Total unique 9-mers (MANE, summed across genes): {mane_kmer_total:>10,}", file=sys.stderr)
    print(f"  Total unique 9-mers (GCF,  summed across genes): {gcf_kmer_total:>10,}", file=sys.stderr)
    print(f"  Done.", file=sys.stderr)


if __name__ == '__main__':
    main()
