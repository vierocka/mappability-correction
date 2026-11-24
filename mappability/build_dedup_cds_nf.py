#!/usr/bin/env python3
"""
Standalone CDS extraction for the nf (Ensembl GRCh38.116) dedup-CDS reference —
run this locally instead of relying on Step 0d inside mappability_correction_
v7_a_nf.sh, which only fires when the whole v7_a_nf pipeline (read simulation +
STAR alignment) runs on the cluster.

Two scopes, matching the two separate dedup-CDS references the _nf pipeline
uses (see Map_scripts/README3.md):
  --scope genome-wide (default) — every transcript_biotype "protein_coding"
      transcript, matching the old v2/v7/v8_b dedup-CDS reference's scope
      exactly (also genome-wide) — for v7_nf/v8_b_nf.
  --scope mane — only tag "MANE_Select" transcripts, a separate, smaller
      reference distinct from both old v2 and genome-wide v7_nf — for v2_nf.

Parses the GTF for CDS features in the chosen scope, splices multi-exon CDS
per transcript via pysam.FastaFile, reverse-complements minus-strand
transcripts, and writes a FASTA with plain versioned-transcript-ID headers
(e.g. ENST00000256078.10 — no '|' characters, unlike NCBI's lcl|... CDS IDs,
which is what let the old v7 pipeline's known ID-grouping bug happen — see
mappability_v7_fix_uf.sh).

If seqkit is on PATH, also runs the dedup step (seqkit rmdup -s) and produces
the final FASTA directly — scp just that one file to the cluster's Ref/
directory (nf reference files live alongside the original NCBI RefSeq ones
there, not a separate nfRef/) and every _nf script's Step 0d will skip
straight past it.
"""

import argparse
import subprocess
import sys
from pathlib import Path
from collections import defaultdict

import pysam

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(RC)[::-1]


def parse_gtf_attrs(col9):
    attrs = defaultdict(list)
    for field in col9.strip().split(';'):
        field = field.strip()
        if not field:
            continue
        key, _, val = field.partition(' ')
        attrs[key].append(val.strip('"'))
    return attrs


def extract_cds(gtf_path, genome_fa_path, raw_out_path, scope):
    mane_select_tx = None
    if scope == 'mane':
        print("Scanning GTF for MANE_Select-tagged transcripts...", flush=True)
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
        print(f"  MANE_Select transcripts: {len(mane_select_tx):,}", flush=True)

    print("Scanning GTF for CDS features...", flush=True)
    tx_cds = defaultdict(list)
    n_cds_lines = 0
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'CDS':
                continue
            a = parse_gtf_attrs(f[8])
            tid = a.get('transcript_id', [None])[0]
            tver = a.get('transcript_version', [None])[0]
            if tid is None:
                continue
            versioned_tid = f"{tid}.{tver}" if tver else tid
            if scope == 'mane':
                if versioned_tid not in mane_select_tx:
                    continue
            else:
                biotype = (a.get('transcript_biotype') or a.get('gene_biotype') or [''])[0]
                if biotype != 'protein_coding':
                    continue
            chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
            tx_cds[versioned_tid].append((chrom, start, end, strand))
            n_cds_lines += 1

    for tid in tx_cds:
        tx_cds[tid].sort(key=lambda e: e[1])

    print(f"  CDS feature lines: {n_cds_lines:,}")
    print(f"  Transcripts with CDS ({scope}): {len(tx_cds):,}", flush=True)

    print("Fetching + splicing CDS sequences...", flush=True)
    fa = pysam.FastaFile(str(genome_fa_path))
    n_written = n_skip = 0
    report_every = 20000
    with open(raw_out_path, 'w') as out:
        for i, (tid, segs) in enumerate(tx_cds.items(), start=1):
            strand = segs[0][3]
            try:
                seq = ''.join(fa.fetch(c, s, e) for (c, s, e, _st) in segs).upper()
            except (KeyError, ValueError):
                n_skip += 1
                continue
            if not seq:
                n_skip += 1
                continue
            if strand == '-':
                seq = revcomp(seq)
            out.write(f">{tid}\n{seq}\n")
            n_written += 1
            if i % report_every == 0:
                print(f"  ...{i:,} / {len(tx_cds):,} transcripts processed", flush=True)
    fa.close()

    print(f"  CDS sequences written (pre-dedup): {n_written:,}  (skipped: {n_skip})")
    print(f"  Saved raw: {raw_out_path}", flush=True)
    return n_written


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gtf', required=True, help='Path to Homo_sapiens.GRCh38.116.gtf')
    ap.add_argument('--genome-fa', required=True, help='Path to Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa')
    ap.add_argument('--scope', choices=['genome-wide', 'mane'], default='genome-wide',
                     help='genome-wide = all protein_coding transcripts (matches old v2/v7/v8_b scope, for v7_nf/v8_b_nf); '
                          'mane = MANE_Select only (separate, smaller reference, for v2_nf)')
    ap.add_argument('--raw-out', default=None, help='Pre-dedup output FASTA (default depends on --scope)')
    ap.add_argument('--dedup-out', default=None, help='Post-dedup output FASTA (default depends on --scope; only written if seqkit is available)')
    ap.add_argument('--threads', type=int, default=4, help='Threads for seqkit rmdup')
    args = ap.parse_args()

    gtf_path = Path(args.gtf)
    genome_fa_path = Path(args.genome_fa)
    suffix = '_mane' if args.scope == 'mane' else ''
    raw_out = Path(args.raw_out or f'Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds{suffix}_raw.fna')
    dedup_out = Path(args.dedup_out or f'Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds{suffix}.fna')

    if not genome_fa_path.with_suffix(genome_fa_path.suffix + '.fai').exists():
        print("Indexing genome with pysam (one-time, builds .fai)...", flush=True)
        pysam.faidx(str(genome_fa_path))

    extract_cds(gtf_path, genome_fa_path, raw_out, args.scope)

    seqkit = subprocess.run(['which', 'seqkit'], capture_output=True, text=True).stdout.strip()
    if not seqkit:
        print("\nseqkit not found on PATH — stopping here.")
        print(f"  Run manually:  seqkit rmdup -s -i {raw_out} -o {dedup_out}")
        print(f"  Then scp {dedup_out} to Ref/{dedup_out.name} on the cluster.")
        sys.exit(0)

    print(f"\nDeduplicating exact-sequence duplicates with seqkit rmdup -s (found at {seqkit})...", flush=True)
    subprocess.run(
        ['seqkit', 'rmdup', '-s', '-i', str(raw_out), '-o', str(dedup_out), '--threads', str(args.threads)],
        check=True,
    )
    print(f"\nDone. scp {dedup_out} to Ref/{dedup_out.name} on the cluster —")
    print("every _nf script's Step 0d will then detect it and skip straight past extraction+dedup.")


if __name__ == '__main__':
    main()
