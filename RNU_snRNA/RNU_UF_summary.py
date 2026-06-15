#!/usr/bin/env python3
"""
Extract UF values for all 23 MANE RNU genes across completed simulation settings.
IDs are parsed from RNU_MANE.fna (ground truth).
Writes RNU_UF_summary.tsv.
"""

import csv
from pathlib import Path

HERE    = Path(__file__).resolve().parent
RESULTS = HERE.parent / "mappability" / "results"

# ── Parse gene → (transcript_id, length, family) from FASTA ──────────────────
def family_of(gene):
    # check longer/specific names before short prefixes to avoid RNU11→RNU1
    for fam in ('RNU11','RNU12','RNU1','RNU2','RNU4','RNU5','RNU6','RNU7'):
        if gene.startswith(fam):
            return fam
    return 'other'

RNU_GENES = {}   # gene_symbol → (enst, length_bp, family)
with open(HERE / "RNU_MANE.fna") as fh:
    for line in fh:
        if line.startswith('>'):
            parts  = line[1:].strip().split('|')
            gene   = parts[0]
            enst   = parts[1]
            length = int(parts[2].replace('bp', ''))
            RNU_GENES[gene] = (enst, length, family_of(gene))

ENST_TO_GENE = {v[0]: k for k, v in RNU_GENES.items()}

# ── Simulation settings to include ───────────────────────────────────────────
SETTINGS = [
    # (column_label, tsv_filename, short_description)
    ('SE_75bp',
     'transcript_uniqueness_factors_MANE_RNA_SE_L75bp.tsv',
     'MANE ref, MANE RNA, SE 75 bp'),
    ('SE_100bp',
     'transcript_uniqueness_factors_MANE_RNA_SE_L100bp.tsv',
     'MANE ref, MANE RNA, SE 100 bp'),
    ('SE_200bp',
     'transcript_uniqueness_factors_MANE_RNA_SE_L200bp.tsv',
     'MANE ref, MANE RNA, SE 200 bp'),
    ('exflank_PE_100bp',
     'transcript_uniqueness_factors_MANE_exonflank_PE_L100bp.tsv',
     'MANE ref, exon-flank reads, PE 100 bp'),
]

# ── Load UF values ────────────────────────────────────────────────────────────
uf_data = {gene: {} for gene in RNU_GENES}

for col_label, fname, _ in SETTINGS:
    path = RESULTS / fname
    if not path.exists():
        for gene in RNU_GENES:
            uf_data[gene][col_label] = 'file_missing'
        continue
    found = set()
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row['transcript_id'] in ENST_TO_GENE:
                gene = ENST_TO_GENE[row['transcript_id']]
                uf_data[gene][col_label] = float(row['uniqueness_factor'])
                found.add(gene)
    for gene in RNU_GENES:
        if col_label not in uf_data[gene]:
            uf_data[gene][col_label] = 'too_short'

# ── Write TSV ─────────────────────────────────────────────────────────────────
col_labels = [s[0] for s in SETTINGS]
OUT = HERE / "RNU_UF_summary.tsv"

with open(OUT, 'w', newline='') as fh:
    writer = csv.writer(fh, delimiter='\t')
    writer.writerow(['gene_symbol', 'family', 'length_bp', 'transcript_id'] + col_labels)
    sort_key = lambda item: (family_of(item[0]), item[0])
    for gene, (enst, length, family) in sorted(RNU_GENES.items(), key=sort_key):
        row = [gene, family, length, enst]
        for col in col_labels:
            v = uf_data[gene].get(col, 'missing')
            row.append(f'{v:.4f}' if isinstance(v, float) else v)
        writer.writerow(row)

print(f"Saved → {OUT}")

# ── Print readable table ──────────────────────────────────────────────────────
hdr_short = [s[0] for s in SETTINGS]
print(f"\n{'Gene':<12} {'Fam':<5} {'Len':>4}  " + "  ".join(f'{h:>16}' for h in hdr_short))
print("-" * (25 + 18 * len(SETTINGS)))
for gene, (enst, length, family) in sorted(RNU_GENES.items(), key=sort_key):
    vals = []
    for col in col_labels:
        v = uf_data[gene].get(col, 'missing')
        vals.append(f'{v:>16.4f}' if isinstance(v, float) else f'{v:>16}')
    print(f"{gene:<12} {family:<5} {length:>4}  " + "  ".join(vals))
