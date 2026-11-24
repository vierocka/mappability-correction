#!/usr/bin/env python3
"""Per-transcript uniqueness factor from a read-back alignment.

funique = (reads from this transcript that re-align uniquely) / (reads simulated from it).

"Uniquely" means MAPQ 255, which is what STAR assigns to a read with exactly one reported
alignment. Read names carry the source transcript (see simulate_reads.py), so the factor is
attributed without needing the alignment to land on any particular feature - the measurement is
of the reference's ability to place the read at all, not of counting.

Only R1 is counted, so a pair contributes once.

usage: uniqueness_from_bam.py --bam F --out F [--mapq-unique 255]
"""
import argparse
import subprocess
from collections import defaultdict


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--mapq-unique", type=int, default=255)
    a = p.parse_args()

    n_pos = defaultdict(int)
    n_uniq = defaultdict(int)
    n_multi = defaultdict(int)
    n_unmapped = defaultdict(int)
    seen = set()

    # -F 256 -F 2048: ignore secondary and supplementary records, count each read once
    cmd = ["samtools", "view", "-F", "0x900", a.bam]
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20) as proc:
        for line in proc.stdout:
            f = line.split("\t", 6)
            name, flag, mapq = f[0], int(f[1]), int(f[4])
            if flag & 0x80:                     # second in pair -> skip, R1 only
                continue
            key = (name, flag & 0x40)
            if key in seen:                     # defensive: duplicate primary record
                continue
            seen.add(key)
            tid = name.rsplit("|", 1)[0]
            n_pos[tid] += 1
            if flag & 0x4:
                n_unmapped[tid] += 1
            elif mapq >= a.mapq_unique:
                n_uniq[tid] += 1
            else:
                n_multi[tid] += 1
        if proc.wait() != 0:
            raise SystemExit("samtools view failed")

    with open(a.out, "w") as fh:
        fh.write("transcript_id\tn_positions\tn_unique\tn_multi\tn_unmapped\tfunique\n")
        for tid in sorted(n_pos):
            tot = n_pos[tid]
            fh.write(f"{tid}\t{tot}\t{n_uniq[tid]}\t{n_multi[tid]}\t{n_unmapped[tid]}\t"
                     f"{n_uniq[tid] / tot:.6f}\n")

    tot = sum(n_pos.values())
    uq = sum(n_uniq.values())
    if not tot:
        raise SystemExit("no reads found in the BAM. Either the alignment is empty, or the read "
                         "names are not of the form <transcript_id>|<position> that "
                         "simulate_reads.py writes.")
    print(f"{len(n_pos)} transcripts, {tot} reads, global funique {uq/tot:.4f}")


if __name__ == "__main__":
    main()
