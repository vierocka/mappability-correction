#!/usr/bin/env python3
"""
Double-check ID conversions between the existing MANE-based pipeline and the
new nf (Ensembl GRCh38.116) reference, before trusting any mapping-difference
result as being caused by the reference rather than an annotation mismatch.

Compares the MANE_Select-tagged transcript set in the nf GTF against
results/id_conversion_table.tsv (built from MANE.GRCh38.v1.5 + NCBI RefSeq).

Does not modify id_conversion_table.tsv or any existing file — writes a new
report at results/id_conversion_check_nf.tsv.
"""

import argparse, csv
from collections import defaultdict
from pathlib import Path


def parse_gtf_attrs(col9):
    """Ensembl GTF attributes: key "value"; repeated keys (esp. 'tag') accumulate."""
    attrs = defaultdict(list)
    for field in col9.strip().split(';'):
        field = field.strip()
        if not field:
            continue
        key, _, val = field.partition(' ')
        attrs[key].append(val.strip('"'))
    return attrs


def load_mane_select_from_gtf(gtf_path):
    """Return {versioned_transcript_id: (versioned_gene_id, gene_name)}."""
    out = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'transcript':
                continue
            a = parse_gtf_attrs(f[8])
            if 'MANE_Select' not in a.get('tag', []):
                continue
            tid  = a.get('transcript_id', [None])[0]
            tver = a.get('transcript_version', [None])[0]
            gid  = a.get('gene_id', [None])[0]
            gver = a.get('gene_version', [None])[0]
            gname = a.get('gene_name', [''])[0]
            if tid is None:
                continue
            versioned_tid = f"{tid}.{tver}" if tver else tid
            versioned_gid = f"{gid}.{gver}" if gid and gver else gid
            out[versioned_tid] = (versioned_gid, gname)
    return out


def load_old_id_table(id_table_path):
    """Return {versioned_transcript_id: (versioned_gene_id, gene_symbol)}."""
    out = {}
    with open(id_table_path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            tid = row['ensembl_transcript_id']
            if not tid:
                continue
            out[tid] = (row['ensembl_gene_id'], row['gene_symbol'])
    return out


def base_id(versioned_id):
    return versioned_id.split('.')[0] if versioned_id else versioned_id


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--gtf', required=True, help='Path to Homo_sapiens.GRCh38.116.gtf (now under Ref/)')
    ap.add_argument('--id-table', default='id_conversion_table.tsv',
                     help='Relative to CWD by default; pass a full path if run from elsewhere')
    ap.add_argument('--out', default='id_conversion_check_nf.tsv',
                     help='Relative to CWD by default; parent dir is created if missing')
    args = ap.parse_args()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading MANE_Select transcripts from nf GTF: {args.gtf}", flush=True)
    new_tx = load_mane_select_from_gtf(args.gtf)
    print(f"  nf MANE_Select transcripts: {len(new_tx):,}")

    print(f"Loading existing id conversion table: {args.id_table}", flush=True)
    old_tx = load_old_id_table(args.id_table)
    print(f"  existing table transcripts: {len(old_tx):,}")

    old_by_base = {base_id(t): t for t in old_tx}
    new_by_base = {base_id(t): t for t in new_tx}

    all_bases = set(old_by_base) | set(new_by_base)
    rows = []
    counts = defaultdict(int)

    for base in sorted(all_bases):
        old_full = old_by_base.get(base)
        new_full = new_by_base.get(base)

        if old_full is None:
            status = 'missing_in_old_table'
            old_gid = old_sym = ''
            new_gid, new_name = new_tx[new_full]
        elif new_full is None:
            status = 'missing_in_nf_gtf'
            new_gid = new_name = ''
            old_gid, old_sym = old_tx[old_full]
        else:
            old_gid, old_sym = old_tx[old_full]
            new_gid, new_name = new_tx[new_full]
            if old_full != new_full:
                status = 'transcript_version_mismatch'
            elif base_id(old_gid) != base_id(new_gid):
                status = 'gene_id_mismatch'
            else:
                status = 'match'

        counts[status] += 1
        rows.append({
            'transcript_id_base': base,
            'old_transcript_id': old_full or '',
            'new_transcript_id': new_full or '',
            'old_ensembl_gene_id': old_gid,
            'new_ensembl_gene_id': new_gid,
            'old_gene_symbol': old_sym if old_full else '',
            'new_gene_name': new_name if new_full else '',
            'status': status,
        })

    with open(args.out, 'w', newline='') as fh:
        fields = ['transcript_id_base', 'old_transcript_id', 'new_transcript_id',
                   'old_ensembl_gene_id', 'new_ensembl_gene_id',
                   'old_gene_symbol', 'new_gene_name', 'status']
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    print("\n── Summary ─────────────────────────────────────────")
    for status, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {status:<28} {n:>7,}")
    print(f"\nWrote: {args.out}")


if __name__ == '__main__':
    main()
