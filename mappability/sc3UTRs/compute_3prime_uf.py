#!/usr/bin/env python3
"""
Compute 3' window uniqueness factors from a simulated-read STAR BAM file.

10x Genomics 3' scRNA-seq captures only the last ~200-600 bp of each
transcript. This script asks: for reads originating from the 3' end only,
what fraction maps uniquely? This is the relevant uniqueness measure for
10x data, not the whole-transcript UF computed by the main pipeline.

Low 3' UF reflects 3' end CONSERVATION — sequences shared between paralogs,
recently duplicated gene family members, or pseudogenes that originated by
retrotransposition of a parental coding 3' end. It is not a quality measure
of the gene itself but a consequence of evolutionary 3' sequence sharing.

Four 3' windows are tested (last N bp of transcript): 250, 400, 500, 600 bp.

Read name format expected: transcript_id|position  (MANE RNA reads)
                       or  transcript_id|exon_idx|position  (exon-flank reads)
Both are handled by rfind('|') to extract the position suffix.

n_positions per transcript is loaded from the existing global UF TSV for
the same run — avoids re-reading the large FASTQ files.

Usage:
    python compute_3prime_uf.py <bam> <ref_tsv> <read_len> <output_tsv> <label>

Arguments:
    bam         STAR BAM (sorted by coordinate, indexed)
    ref_tsv     Existing per-transcript UF TSV for this run (provides n_positions)
    read_len    Read length used in this run (bp)
    output_tsv  Output path
    label       Run label, e.g. v1_c_genome_SE_100bp
"""

import sys
import csv
import pysam
from collections import defaultdict
from pathlib import Path

WINDOWS = [250, 400, 500, 600]


def main():
    if len(sys.argv) != 6:
        print(__doc__)
        sys.exit(1)

    bam_path = Path(sys.argv[1])
    ref_tsv  = Path(sys.argv[2])
    read_len = int(sys.argv[3])
    out_path = Path(sys.argv[4])
    label    = sys.argv[5]

    print(f"\n{'═'*60}")
    print(f"3' UF analysis: {label}")
    print(f"  BAM     : {bam_path}")
    print(f"  ref TSV : {ref_tsv}")
    print(f"  read_len: {read_len} bp")
    print(f"  windows : {WINDOWS}")
    print(f"{'═'*60}", flush=True)

    # ── Load n_positions and overall UF from existing TSV ─────────────────────
    n_pos      = {}   # {transcript_id: n_positions}
    uf_overall = {}   # {transcript_id: uniqueness_factor}
    with open(ref_tsv) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            tid = row['transcript_id']
            n_pos[tid]      = int(row['n_positions'])
            uf_overall[tid] = float(row['uniqueness_factor'])
    print(f"  Loaded {len(n_pos):,} transcripts from {ref_tsv.name}", flush=True)

    # ── Gene symbols ──────────────────────────────────────────────────────────
    enst_to_gene = {}
    id_table = ref_tsv.parent / 'id_conversion_table.tsv'
    if id_table.exists():
        with open(id_table) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                enst_to_gene[row['ensembl_transcript_id']] = row['gene_symbol']

    # ── Single BAM pass ───────────────────────────────────────────────────────
    # For each R1 read with MAPQ=255, record whether it falls in each 3' window.
    # 3' window W: positions pos >= max_pos - W + 1  (max_pos = n_positions - 1)

    unique_total = defaultdict(int)
    unique_3p    = {W: defaultdict(int) for W in WINDOWS}

    bam = pysam.AlignmentFile(str(bam_path), 'rb')
    n_proc = n_uniq = 0
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.is_read2:
            continue   # count R1 only; consistent with SE comparison

        qname = read.query_name
        pipe  = qname.rfind('|')
        if pipe < 0:
            continue
        try:
            pos = int(qname[pipe + 1:])
        except ValueError:
            continue
        tid = qname[:pipe]

        if tid not in n_pos:
            continue
        n_proc += 1

        if read.mapping_quality == 255:
            n_uniq += 1
            unique_total[tid] += 1
            max_p = n_pos[tid] - 1
            for W in WINDOWS:
                if pos >= max_p - W + 1:
                    unique_3p[W][tid] += 1

        if n_proc % 5_000_000 == 0:
            print(f"  {n_proc // 1_000_000}M reads ...", flush=True)

    bam.close()
    print(f"  R1 reads processed: {n_proc:,}  MAPQ=255: {n_uniq:,}", flush=True)

    # ── Build output rows ─────────────────────────────────────────────────────
    fields = ['transcript_id', 'gene_symbol', 'transcript_len', 'read_len',
              'label', 'n_positions', 'uf_overall']
    for W in WINDOWS:
        fields += [f'n_sim_3p{W}', f'n_unique_3p{W}', f'uf_3p{W}', f'delta_3p{W}']

    rows = []
    for tid in sorted(n_pos):
        np_val  = n_pos[tid]
        tx_len  = np_val + read_len - 1
        uf_ov   = uf_overall[tid]
        row = {
            'transcript_id':  tid,
            'gene_symbol':    enst_to_gene.get(tid, ''),
            'transcript_len': tx_len,
            'read_len':       read_len,
            'label':          label,
            'n_positions':    np_val,
            'uf_overall':     round(uf_ov, 6),
        }
        for W in WINDOWS:
            # Denominator: reads simulated in the 3' window.
            # For transcripts shorter than W, the entire transcript is the window.
            n_sim   = min(W, np_val)
            n_uniq_w = unique_3p[W].get(tid, 0)
            uf_3p   = round(n_uniq_w / n_sim, 6) if n_sim > 0 else 0.0
            row[f'n_sim_3p{W}']    = n_sim
            row[f'n_unique_3p{W}'] = n_uniq_w
            row[f'uf_3p{W}']       = uf_3p
            row[f'delta_3p{W}']    = round(uf_3p - uf_ov, 6)
        rows.append(row)

    # Sort by 500 bp window delta (most 3'-degraded first)
    rows.sort(key=lambda r: r['delta_3p500'])

    with open(out_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    # ── Summary stats ─────────────────────────────────────────────────────────
    n = len(rows)
    for W in WINDOWS:
        ufs    = [r[f'uf_3p{W}']    for r in rows]
        deltas = [r[f'delta_3p{W}'] for r in rows]
        n_worse = sum(1 for d in deltas if d < -0.10)
        n_zero  = sum(1 for f in ufs  if f == 0.0)
        print(f"  [{W:3d} bp]  mean_uf={sum(ufs)/n:.4f}  "
              f"uf=0: {n_zero:5,}  delta<-0.10: {n_worse:5,}")
    print(f"\n  Saved → {out_path}", flush=True)


if __name__ == '__main__':
    main()
