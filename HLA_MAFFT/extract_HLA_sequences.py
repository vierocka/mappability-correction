#!/usr/bin/env python3
"""
Extract HLA class I sequences from three sources:
  1. MANE RNA FASTA       → HLA_MANE.fna        (spliced mRNA, MANE Select transcripts)
  2. Dedup CDS FASTA      → HLA_dedup_cds.fna   (CDS sequences, deduplicated alleles)
  3. Genome + GFF         → HLA_genomic.fna      (gene body incl. introns + pseudogenes)

Genes targeted:
  Protein-coding : HLA-A, HLA-B, HLA-C, HLA-E, HLA-F, HLA-G
  Pseudogenes    : HLA-V, HLA-P, HLA-H, HLA-T, HLA-K, HLA-U,
                   HLA-W, HLA-J, HLA-L, HLA-N, MICE, MICG, MICF, MICD
"""

import subprocess, re
from pathlib import Path

REF = Path("/home/veve/Dropbox/Self-Nonself/Reference/mappability-correction/Ref")
OUT = Path("/home/veve/Dropbox/Self-Nonself/Reference/mappability-correction/HLA_MAFFT")
OUT.mkdir(parents=True, exist_ok=True)

HLA_CODING    = {'HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G'}
HLA_PSEUDO    = {'HLA-V', 'HLA-P', 'HLA-H', 'HLA-T', 'HLA-K', 'HLA-U',
                 'HLA-W', 'HLA-J', 'HLA-L', 'HLA-N', 'MICE', 'MICG', 'MICF', 'MICD'}
ALL_HLA       = HLA_CODING | HLA_PSEUDO

# ── Version 1: MANE RNA ───────────────────────────────────────────────────────
print("=" * 60)
print("Version 1: MANE spliced mRNA")
print("=" * 60)
mane_rna = REF / "MANE.GRCh38.v1.5.ensembl_rna.fna"
out_mane  = OUT / "HLA_MANE.fna"
n = 0
with open(mane_rna) as fh, open(out_mane, 'w') as out:
    write = False
    for line in fh:
        if line.startswith('>'):
            m = re.search(r'gene_symbol:(\S+)', line)
            gene = m.group(1) if m else ''
            write = gene in ALL_HLA
            if write:
                # Rename header: >GENE_NAME original_header
                out.write(f">{gene} {line[1:].strip()}\n")
                print(f"  {gene}")
                n += 1
                continue
        if write:
            out.write(line)
print(f"  Total sequences: {n}  →  {out_mane.name}\n")

# ── Version 2: Deduplicated CDS ───────────────────────────────────────────────
print("=" * 60)
print("Version 2: Deduplicated CDS")
print("=" * 60)
dedup_cds = REF / "GCF_000001405.40_GRCh38.p14_cds_from_genomic.dedup.fna"
out_dedup  = OUT / "HLA_dedup_cds.fna"
n = 0
counts = {}
with open(dedup_cds) as fh, open(out_dedup, 'w') as out:
    write = False
    gene_cur = ''
    for line in fh:
        if line.startswith('>'):
            m = re.search(r'\[gene=([^\]]+)\]', line)
            gene_cur = m.group(1) if m else ''
            write = gene_cur in ALL_HLA
            if write:
                counts[gene_cur] = counts.get(gene_cur, 0) + 1
                idx = counts[gene_cur]
                # Rename header for readability
                pid_m = re.search(r'\[protein_id=([^\]]+)\]', line)
                pid   = pid_m.group(1) if pid_m else f"seq{idx}"
                out.write(f">{gene_cur}_{idx} {pid}\n")
                continue
        if write:
            out.write(line)
for gene, cnt in sorted(counts.items()):
    print(f"  {gene}: {cnt} sequences")
    n += cnt
print(f"  Total sequences: {n}  →  {out_dedup.name}\n")

# ── Version 3: Genomic gene body (GFF + samtools faidx) ──────────────────────
print("=" * 60)
print("Version 3: Genomic gene body (coding + pseudogenes)")
print("=" * 60)
gff_path  = REF / "GCF_000001405.40_GRCh38.p14_genomic.gff"
genome_fa = REF / "GCF_000001405.40_GRCh38.p14_genomic.fna"
out_geno  = OUT / "HLA_primaryChr.fna"

# Parse GFF — keep one entry per gene (longest if duplicates on alt contigs)
gene_coords = {}
with open(gff_path) as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('gene', 'pseudogene'): continue
        m = re.search(r'Name=([^;]+)', f[8])
        if not m: continue
        gene = m.group(1)
        if gene not in ALL_HLA: continue
        chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
        # Prefer primary chromosome; otherwise take longest entry
        cur = gene_coords.get(gene)
        primary = chrom == 'NC_000006.12'
        cur_primary = cur is not None and cur[0] == 'NC_000006.12'
        if cur is None or (primary and not cur_primary) or \
           (primary == cur_primary and (end - start) > (cur[2] - cur[1])):
            gene_coords[gene] = (chrom, start, end, strand)

print(f"  Genes found in GFF: {sorted(gene_coords)}\n")
n = 0
with open(out_geno, 'w') as out:
    for gene in sorted(gene_coords):
        chrom, start, end, strand = gene_coords[gene]
        region = f"{chrom}:{start + 1}-{end}"
        try:
            res = subprocess.run(
                ['samtools', 'faidx', str(genome_fa), region],
                capture_output=True, text=True, check=True)
            lines = res.stdout.strip().split('\n')
            category = 'pseudo' if gene in HLA_PSEUDO else 'coding'
            out.write(f">{gene} {region} strand={strand} [{category}]\n")
            out.write('\n'.join(lines[1:]) + '\n')
            bp = end - start
            print(f"  {gene:10s}  {chrom}:{start+1}-{end}  strand={strand}  {bp:,} bp  [{category}]")
            n += 1
        except subprocess.CalledProcessError:
            print(f"  WARNING: failed to fetch {gene} ({region})")
print(f"\n  Total sequences: {n}  →  {out_geno.name}  (NC_000006.12 primary chr only, one entry per gene)")
