#!/usr/bin/env python3
"""
Generate a table of unique 9-mer peptides per MANE transcript.

Since MANE defines one canonical transcript per gene, each ENST has
exactly one protein (ENSP). The script maps ENST → ENSP from the
protein FASTA headers, then slides a 9-AA window per protein.

Input:
  1. MANE.GRCh38.v1.5.ensembl_protein.faa  (ENSP FASTA with ENST in header)

Output:
  mane_9mer_table.tsv  (gene_id [ENST], peptide_9mer, count_in_gene)

Usage:
    python3 generate_mane_9mer_table.py \
        MANE.GRCh38.v1.5.ensembl_protein.faa \
        mane_9mer_table.tsv
"""

import sys
import re
from collections import defaultdict

fasta_file = sys.argv[1]
out_file = sys.argv[2]

KMER = 9

# ── Parse FASTA: extract ENST → protein sequence ──
# Header format:
# >ENSP00000263100.2 pep ... transcript:ENST00000263100.8 ... gene_symbol:A1BG ...
transcript_re = re.compile(r'transcript:(ENST\d+\.\d+(?:_[A-Za-z0-9]+)?)')   # PIPELINE FIX: keep HLA allele suffix (_A2301)
symbol_re = re.compile(r'gene_symbol:(\S+)')

transcripts = {}  # ENST -> { 'seq': str, 'symbol': str }
current_enst = None
current_symbol = None
current_seq = []

with open(fasta_file) as f:
    for line in f:
        if line.startswith('>'):
            # Save previous
            if current_enst:
                transcripts[current_enst] = {
                    'seq': ''.join(current_seq).rstrip('*'),
                    'symbol': current_symbol or 'NA',
                }
            # Parse header
            tm = transcript_re.search(line)
            sm = symbol_re.search(line)
            current_enst = tm.group(1) if tm else None
            current_symbol = sm.group(1) if sm else None
            current_seq = []
        else:
            current_seq.append(line.strip())

    # Last entry
    if current_enst:
        transcripts[current_enst] = {
            'seq': ''.join(current_seq).rstrip('*'),
            'symbol': current_symbol or 'NA',
        }

print(f"Parsed {len(transcripts)} MANE transcripts from FASTA")

# ── Generate 9-mers per transcript ──
total_rows = 0

with open(out_file, 'w') as out:
    out.write("gene_id\tgene_symbol\tpeptide_9mer\tcount_in_gene\n")

    for enst in sorted(transcripts.keys()):
        seq = transcripts[enst]['seq']
        symbol = transcripts[enst]['symbol']

        kmer_counts = defaultdict(int)
        for i in range(len(seq) - KMER + 1):
            kmer = seq[i:i + KMER]
            if 'X' in kmer or '*' in kmer or 'U' in kmer:
                continue
            kmer_counts[kmer] += 1

        for kmer in sorted(kmer_counts.keys()):
            out.write(f"{enst}\t{symbol}\t{kmer}\t{kmer_counts[kmer]}\n")
            total_rows += 1

print(f"Wrote {total_rows} transcript-peptide entries to {out_file}")
print(f"Covering {len(transcripts)} transcripts")
