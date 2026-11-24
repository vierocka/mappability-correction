#!/usr/bin/env python3
"""Tile transcripts into error-free reads for a mappability read-back test.

For every transcript, a read is placed at each position with a given step. With --mode PE a
mate is emitted as the reverse complement of the fragment's far end, so the pair spans a
realistic insert rather than the same bases twice.

Read names encode the source transcript, so uniqueness can be attributed after alignment:
    <transcript_id>|<position>

usage: simulate_reads.py --fasta F --r1 OUT.fq.gz [--r2 OUT.fq.gz] [--length 100] [--step 1]
                         [--mode SE|PE] [--insert 300] [--flank-fasta F]
"""
import argparse
import gzip

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
rc = lambda s: s.translate(COMP)[::-1]


def read_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--r1", required=True)
    p.add_argument("--r2")
    p.add_argument("--length", type=int, default=100)
    p.add_argument("--step", type=int, default=1)
    p.add_argument("--mode", choices=["SE", "PE"], default="SE")
    p.add_argument("--insert", type=int, default=300)
    a = p.parse_args()

    if a.mode == "PE" and not a.r2:
        raise SystemExit("--mode PE needs --r2")

    qual = "I" * a.length
    n_reads = n_tx = n_short = 0
    f1 = gzip.open(a.r1, "wt")
    f2 = gzip.open(a.r2, "wt") if a.mode == "PE" else None
    try:
        for tid, seq in read_fasta(a.fasta):
            seq = seq.upper()
            if len(seq) < a.length:
                n_short += 1
                continue
            n_tx += 1
            # fragment end is bounded by the transcript; short transcripts get a shorter insert
            for i in range(0, len(seq) - a.length + 1, a.step):
                r1 = seq[i:i + a.length]
                name = f"{tid}|{i}"
                f1.write(f"@{name}/1\n{r1}\n+\n{qual}\n")
                if f2 is not None:
                    frag_end = min(i + a.insert, len(seq))
                    start2 = max(frag_end - a.length, i)
                    f2.write(f"@{name}/2\n{rc(seq[start2:frag_end]).ljust(a.length, 'N')[:a.length]}\n+\n{qual}\n")
                n_reads += 1
    finally:
        f1.close()
        if f2 is not None:
            f2.close()

    print(f"{n_reads} reads ({a.mode} {a.length} bp, step {a.step}) from {n_tx} transcripts; "
          f"{n_short} shorter than the read length were skipped")
    if n_reads == 0:
        raise SystemExit(f"no reads simulated: {a.fasta} held no sequence at least "
                         f"{a.length} bp long. An empty or truncated transcript FASTA is the "
                         f"usual cause - check the rule that built it.")


if __name__ == "__main__":
    main()
