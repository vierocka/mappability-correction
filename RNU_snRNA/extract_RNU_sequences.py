#!/usr/bin/env python3
"""
Extract all RNU* sequences from the MANE Select v1.5 RNA FASTA.

Outputs:
  RNU_MANE.fna         – all 23 sequences, header: >gene_symbol|transcript_id|lengthbp
  RNU_MANE_unique.fna  – same sequences with identical copies collapsed (seqkit rmdup)
  per-family FASTAs    – RNU1_family.fna, RNU2_family.fna, RNU4_family.fna,
                         RNU5_family.fna, RNU6_family.fna
                         (families with ≥2 members; RNU7-1, RNU11, RNU12 are singletons)
"""

from pathlib import Path

HERE  = Path(__file__).resolve().parent
MANE  = HERE.parent / "Ref" / "MANE.GRCh38.v1.5.ensembl_rna.fna"
OUT   = HERE / "RNU_MANE.fna"

FAMILIES = {
    'RNU1': ['RNU1-1', 'RNU1-2', 'RNU1-3', 'RNU1-4'],
    'RNU2': ['RNU2-1', 'RNU2-2', 'RNU2-63P'],
    'RNU4': ['RNU4-1', 'RNU4-2', 'RNU4ATAC'],
    'RNU5': ['RNU5A-1', 'RNU5B-1', 'RNU5D-1', 'RNU5E-1'],
    'RNU6': ['RNU6-1', 'RNU6-2', 'RNU6-7', 'RNU6-8', 'RNU6-9', 'RNU6ATAC'],
}

seqs = []
current_header = None
current_seq    = []

with open(MANE) as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith('>'):
            if current_header and 'gene_symbol:RNU' in current_header:
                seq = ''.join(current_seq)
                parts = {kv.split(':')[0]: ':'.join(kv.split(':')[1:])
                         for kv in current_header[1:].split() if ':' in kv}
                gene = parts.get('gene_symbol', '?')
                tid  = current_header[1:].split()[0]
                seqs.append((gene, tid, seq))
            current_header = line
            current_seq    = []
        else:
            current_seq.append(line)

if current_header and 'gene_symbol:RNU' in current_header:
    seq   = ''.join(current_seq)
    parts = {kv.split(':')[0]: ':'.join(kv.split(':')[1:])
             for kv in current_header[1:].split() if ':' in kv}
    gene  = parts.get('gene_symbol', '?')
    tid   = current_header[1:].split()[0]
    seqs.append((gene, tid, seq))

seqs.sort(key=lambda x: x[0])

# ── write RNU_MANE.fna ────────────────────────────────────────────────────────
with open(OUT, 'w') as fh:
    for gene, tid, seq in seqs:
        fh.write(f'>{gene}|{tid}|{len(seq)}bp\n{seq}\n')

print(f"Extracted {len(seqs)} RNU sequences → {OUT}")
for gene, tid, seq in seqs:
    print(f"  {gene:12s}  {tid}  {len(seq)} bp")

# ── write per-family FASTAs ───────────────────────────────────────────────────
seq_dict = {gene: (tid, seq) for gene, tid, seq in seqs}

for fam, members in FAMILIES.items():
    out_fam = HERE / f"{fam}_family.fna"
    with open(out_fam, 'w') as fh:
        for gene in members:
            if gene in seq_dict:
                tid, seq = seq_dict[gene]
                fh.write(f'>{gene}|{tid}|{len(seq)}bp\n{seq}\n')
    present = [g for g in members if g in seq_dict]
    print(f"\n{fam}_family.fna → {len(present)} sequences: {', '.join(present)}")
