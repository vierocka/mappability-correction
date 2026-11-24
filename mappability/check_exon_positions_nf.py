#!/usr/bin/env python3
"""
Double-check exon positions between the existing MANE GFF (matched to the NCBI
RefSeq genome) and the new nf Ensembl GTF (MANE_Select-tagged transcripts),
before trusting mapping differences as reference-driven rather than an
annotation/patch-level coordinate shift.

Both files describe GRCh38, but the NCBI RefSeq primary assembly (p14)
includes reference patches that Ensembl's plain "primary_assembly" FASTA does
not — coordinates for a given transcript should match after chromosome-name
translation, but patched regions are exactly where they might not.

Only needs the nf GTF (not the nf genome FASTA or STAR) — independent of
mappability_correction_v1_a_nf.sh. Chromosome-name translation (old "chr1"
style -> new "1" style) is derived on the fly from the existing
chr_name_map.tsv using the same canonical rule that script uses; if an
expanded chr_name_map_nf.tsv already exists (because that script has since
been run), pass --chr-map to use its scaffold-aware version instead.

All paths default to plain filenames resolved relative to the current working
directory (matching how this is invoked on the cluster) — pass explicit paths
if running from elsewhere. Writes a new report (default
exon_position_check_nf.tsv) — does not modify any existing file.
"""

import argparse, csv
from collections import defaultdict
from pathlib import Path


def parse_gff3_attrs(col9):
    attrs = {}
    for field in col9.split(';'):
        if '=' in field:
            k, v = field.split('=', 1)
            attrs[k.strip()] = v.strip()
    return attrs


def parse_gtf_attrs(col9):
    attrs = defaultdict(list)
    for field in col9.strip().split(';'):
        field = field.strip()
        if not field:
            continue
        key, _, val = field.partition(' ')
        attrs[key].append(val.strip('"'))
    return attrs


def canonical_ensembl_name(chr_name):
    """Same rule as Step 0 of mappability_correction_v1_a_nf.sh."""
    if chr_name in ("chrM", "chrMT"):
        return "MT"
    if chr_name.startswith("chr"):
        return chr_name[3:]
    return chr_name


def load_chr_map(path):
    """chr_name -> ensembl_name.

    If `path` is an expanded chr_name_map_nf.tsv (has a 'note' column), use its
    'canonical' rows. Otherwise treat it as a plain chr_name/ncbi_name map (the
    existing results/chr_name_map.tsv) and derive the ensembl name on the fly —
    this is what lets this check run without mappability_correction_v1_a_nf.sh
    having been run first.
    """
    m = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if reader.fieldnames and 'note' in reader.fieldnames:
            for row in reader:
                if row.get('note') == 'canonical':
                    m[row['chr_name']] = row['ensembl_name']
        else:
            fh.seek(0)
            for line in fh:
                parts = line.strip().split('\t')
                if len(parts) != 2:
                    continue
                chr_name, _ncbi_name = parts
                m[chr_name] = canonical_ensembl_name(chr_name)
    return m


def load_old_gff_exons(gff_path):
    """{transcript_id: [(chr_name, start0, end, strand), ...]} sorted by start."""
    tx_exons = defaultdict(list)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'exon':
                continue
            a = parse_gff3_attrs(f[8])
            tid = a.get('transcript_id')
            if tid is None:
                continue
            start = int(f[3]) - 1
            end = int(f[4])
            tx_exons[tid].append((f[0], start, end, f[6]))
    for tid in tx_exons:
        tx_exons[tid].sort(key=lambda e: e[1])
    return tx_exons


def load_new_gtf_exons(gtf_path):
    """{transcript_id: [(chrom, start0, end, strand), ...]} for MANE_Select transcripts only."""
    mane_select_tx = set()
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
            tid = a.get('transcript_id', [None])[0]
            tver = a.get('transcript_version', [None])[0]
            if tid is None:
                continue
            mane_select_tx.add(f"{tid}.{tver}" if tver else tid)

    tx_exons = defaultdict(list)
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'exon':
                continue
            a = parse_gtf_attrs(f[8])
            tid = a.get('transcript_id', [None])[0]
            tver = a.get('transcript_version', [None])[0]
            if tid is None:
                continue
            versioned_tid = f"{tid}.{tver}" if tver else tid
            if versioned_tid not in mane_select_tx:
                continue
            start = int(f[3]) - 1
            end = int(f[4])
            tx_exons[versioned_tid].append((f[0], start, end, f[6]))
    for tid in tx_exons:
        tx_exons[tid].sort(key=lambda e: e[1])
    return tx_exons


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--old-gff', default='MANE.GRCh38.v1.5.ensembl_genomic.gff',
                     help='Relative to CWD by default; pass a full path if run from elsewhere')
    ap.add_argument('--new-gtf', required=True, help='Path to Homo_sapiens.GRCh38.116.gtf (now under Ref/)')
    ap.add_argument('--chr-map', default='chr_name_map.tsv',
                     help='Plain chr_name/ncbi_name map (default) or an expanded '
                          'chr_name_map_nf.tsv if mappability_correction_v1_a_nf.sh has been run')
    ap.add_argument('--out', default='exon_position_check_nf.tsv',
                     help='Relative to CWD by default; parent dir is created if missing')
    args = ap.parse_args()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading chr_name -> ensembl_name map: {args.chr_map}", flush=True)
    chr_map = load_chr_map(args.chr_map)
    print(f"  {len(chr_map)} canonical chromosomes")

    print(f"Loading old MANE GFF exons: {args.old_gff}", flush=True)
    old_exons = load_old_gff_exons(args.old_gff)
    print(f"  {len(old_exons):,} transcripts")

    print(f"Loading new nf GTF exons (MANE_Select only): {args.new_gtf}", flush=True)
    new_exons = load_new_gtf_exons(args.new_gtf)
    print(f"  {len(new_exons):,} transcripts")

    def base_id(t):
        return t.split('.')[0] if t else t

    old_by_base = {base_id(t): t for t in old_exons}
    new_by_base = {base_id(t): t for t in new_exons}
    all_bases = sorted(set(old_by_base) | set(new_by_base))

    rows = []
    counts = defaultdict(int)

    for base in all_bases:
        old_full = old_by_base.get(base)
        new_full = new_by_base.get(base)

        if old_full is None:
            status = 'missing_in_old_gff'
        elif new_full is None:
            status = 'missing_in_new_gtf'
        else:
            old_ex = old_exons[old_full]
            new_ex = new_exons[new_full]
            if len(old_ex) != len(new_ex):
                status = 'exon_count_mismatch'
            else:
                status = 'match'
                mismatch_detail = ''
                for (o_chr, o_start, o_end, o_strand), (n_chr, n_start, n_end, n_strand) in zip(old_ex, new_ex):
                    expected_chr = chr_map.get(o_chr, o_chr)
                    if expected_chr != n_chr:
                        status = 'chrom_mismatch'
                        mismatch_detail = f"old={o_chr}->{expected_chr} new={n_chr}"
                        break
                    if o_start != n_start or o_end != n_end:
                        status = 'coordinate_mismatch'
                        mismatch_detail = f"old={o_chr}:{o_start}-{o_end} new={n_chr}:{n_start}-{n_end}"
                        break
                    if o_strand != n_strand:
                        status = 'strand_mismatch'
                        mismatch_detail = f"old={o_strand} new={n_strand}"
                        break

        counts[status] += 1
        rows.append({
            'transcript_id_base': base,
            'old_transcript_id': old_full or '',
            'new_transcript_id': new_full or '',
            'n_exons_old': len(old_exons.get(old_full, [])) if old_full else '',
            'n_exons_new': len(new_exons.get(new_full, [])) if new_full else '',
            'status': status,
            'detail': mismatch_detail if status not in ('match', 'missing_in_old_gff', 'missing_in_new_gtf', 'exon_count_mismatch') else '',
        })

    with open(args.out, 'w', newline='') as fh:
        fields = ['transcript_id_base', 'old_transcript_id', 'new_transcript_id',
                   'n_exons_old', 'n_exons_new', 'status', 'detail']
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    print("\n── Summary ─────────────────────────────────────────")
    for status, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {status:<24} {n:>7,}")
    print(f"\nWrote: {args.out}")


if __name__ == '__main__':
    main()
