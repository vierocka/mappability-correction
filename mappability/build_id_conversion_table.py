#!/usr/bin/env python3
"""
Build a MANE ↔ RefSeq ID conversion table.

Sources:
  MANE.GRCh38.v1.5.ensembl_genomic.gff  → Ensembl gene/transcript/protein IDs
                                            + RefSeq transcript ID (NM_) via Dbxref
  GCF_000001405.40_GRCh38.p14_genomic.gff → RefSeq GeneID and protein ID (NP_)
                                              linked via NM_ transcript ID

Output columns:
  gene_symbol, ensembl_gene_id, ensembl_transcript_id, ensembl_protein_id,
  refseq_transcript_id, refseq_gene_id, refseq_protein_id
"""

import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
REF  = HERE.parent / "Ref"
OUT  = HERE / "results" / "id_conversion_table.tsv"

MANE_GFF   = REF / "MANE.GRCh38.v1.5.ensembl_genomic.gff"
REFSEQ_GFF = REF / "GCF_000001405.40_GRCh38.p14_genomic.gff"


def parse_attrs(col9):
    attrs = {}
    for field in col9.split(';'):
        if '=' in field:
            k, v = field.split('=', 1)
            attrs[k.strip()] = v.strip()
    return attrs


# ── Step 1: parse MANE GFF — transcript features ──────────────────────────────
print("Parsing MANE GFF ...", flush=True)
mane_rows = {}   # NM_ → {ensembl fields}

with open(MANE_GFF) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'transcript':
            continue
        a = parse_attrs(f[8])
        nm_m = re.search(r'RefSeq:(NM_\S+)', a.get('Dbxref', ''))
        if not nm_m:
            continue
        nm = nm_m.group(1)
        mane_rows[nm] = {
            'gene_symbol':           a.get('gene_name', ''),
            'ensembl_gene_id':       a.get('gene_id', ''),
            'ensembl_transcript_id': a.get('transcript_id', ''),
            'ensembl_protein_id':    a.get('protein_id', ''),
            'refseq_transcript_id':  nm,
            'refseq_gene_id':        '',
            'refseq_protein_id':     '',
        }

print(f"  MANE transcripts with NM_ cross-ref: {len(mane_rows)}")

# ── Step 2: parse RefSeq GFF — mRNA features → GeneID ─────────────────────────
# and CDS features → protein_id (NP_)
print("Parsing RefSeq GFF (this may take a moment) ...", flush=True)

nm_to_geneid  = {}   # NM_ → GeneID string
nm_to_protein = {}   # NM_ → NP_ (first hit per NM_)

with open(REFSEQ_GFF) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        ftype = f[2]

        if ftype == 'mRNA':
            a  = parse_attrs(f[8])
            nm = a.get('transcript_id', '')
            if nm not in mane_rows:
                continue
            gid_m = re.search(r'GeneID:(\d+)', a.get('Dbxref', ''))
            if gid_m:
                nm_to_geneid[nm] = gid_m.group(1)

        elif ftype == 'CDS':
            a      = parse_attrs(f[8])
            parent = a.get('Parent', '')          # rna-NM_XXXXXX.X
            nm_m   = re.match(r'rna-(NM_\S+)', parent)
            if not nm_m:
                continue
            nm = nm_m.group(1)
            if nm not in mane_rows or nm in nm_to_protein:
                continue
            np_ = a.get('protein_id', '')
            if np_:
                nm_to_protein[nm] = np_

print(f"  GeneIDs matched:   {len(nm_to_geneid)}")
print(f"  Protein IDs matched: {len(nm_to_protein)}")

# ── Step 3: join and write ─────────────────────────────────────────────────────
fields = ['gene_symbol', 'ensembl_gene_id', 'ensembl_transcript_id',
          'ensembl_protein_id', 'refseq_transcript_id',
          'refseq_gene_id', 'refseq_protein_id']

n_written = n_missing_gene = n_missing_prot = 0

with open(OUT, 'w') as fh:
    fh.write('\t'.join(fields) + '\n')
    for nm, row in sorted(mane_rows.items(), key=lambda x: x[1]['gene_symbol']):
        row['refseq_gene_id']    = nm_to_geneid.get(nm, '')
        row['refseq_protein_id'] = nm_to_protein.get(nm, '')
        fh.write('\t'.join(row[f] for f in fields) + '\n')
        n_written += 1
        if not row['refseq_gene_id']:    n_missing_gene += 1
        if not row['refseq_protein_id']: n_missing_prot += 1

print(f"\nRows written      : {n_written}")
print(f"Missing GeneID    : {n_missing_gene}")
print(f"Missing protein ID: {n_missing_prot}")
print(f"\nSaved → {OUT}")
