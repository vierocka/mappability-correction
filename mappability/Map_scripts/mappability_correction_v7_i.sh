#!/bin/bash
# v7_i — dedup CDS reads | full genome | SE | 200 bp

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$SCRIPT_DIR/../../Ref"
DEDUP_FA="$REF/GCF_000001405.40_GRCh38.p14_cds_from_genomic.dedup.fna"
GENOME_FA="$REF/GCF_000001405.40_GRCh38.p14_genomic.fna"
GENOME_IDX="$REF/STAR_index_genome"
OUTDIR="$SCRIPT_DIR/../results/mappability_dedup_cds_reads_SE_L200bp"
THREADS=64
READ_LEN=200

mkdir -p "$OUTDIR/simreads" "$OUTDIR/bam" "$OUTDIR/logs"

echo "════════════════════════════════════════════════════════"
echo "Mappability — v7_i (dedup CDS reads, full genome, SE, 200 bp)"
echo "READ_LEN=$READ_LEN  THREADS=$THREADS"
echo "════════════════════════════════════════════════════════"

echo ""
echo "Step 0: Checking STAR index (full genome)..."
mkdir -p "$GENOME_IDX"
if [ ! -f "$GENOME_IDX/SA" ]; then
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$GENOME_IDX" \
        --genomeFastaFiles    "$GENOME_FA" \
        --genomeSAindexNbases 14 \
        --genomeChrBinNbits   14 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR/logs/star_index.log"
else
    echo "Step 0: index already exists, skipping."
fi

echo ""
echo "Step 1: Simulating reads from dedup CDS FASTA (SE, READ_LEN=$READ_LEN)..."

if [ -f "$OUTDIR/simreads/sim_R1.fastq" ]; then
    echo "  Reads exist, skipping simulation."
else

python3 - << PYEOF
from pathlib import Path

cds_fa   = Path("$DEDUP_FA")
outdir   = Path("$OUTDIR")
read_len = int("$READ_LEN")

def parse_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name: yield name, ''.join(seq)
                name = line[1:].split()[0]; seq = []
            else: seq.append(line.upper())
    if name: yield name, ''.join(seq)

r1_path = outdir / 'simreads/sim_R1.fastq'
qual    = 'I' * read_len
n_reads = n_seq = n_skip = 0

with open(r1_path, 'w') as r1:
    for seq_id, seq in parse_fasta(cds_fa):
        L = len(seq)
        if L < read_len: n_skip += 1; continue
        n_seq += 1
        for pos in range(L - read_len + 1):
            fwd = seq[pos : pos + read_len]
            if fwd.count('N') > read_len // 5: continue
            r1.write(f'@{seq_id}|{pos}\n{fwd}\n+\n{qual}\n')
            n_reads += 1

print(f"  CDS sequences: {n_seq:,}  (too short: {n_skip})")
print(f"  Total reads simulated: {n_reads:,}", flush=True)
PYEOF

fi
echo "Step 1 complete."

echo ""
echo "Step 2: STAR alignment (SE, full genome, $THREADS threads)..."

STAR \
    --runThreadN            $THREADS \
    --genomeDir             "$GENOME_IDX" \
    --readFilesIn           "$OUTDIR/simreads/sim_R1.fastq" \
    --outSAMtype            BAM SortedByCoordinate \
    --outSAMattributes      NH HI AS NM \
    --outSAMmultNmax        1 \
    --outFilterMultimapNmax 40 \
    --outBAMsortingThreadN  $THREADS \
    --outBAMsortingBinsN    20 \
    --limitBAMsortRAM       160000000000 \
    --outFileNamePrefix     "$OUTDIR/bam/sim_" \
    2>&1 | tee "$OUTDIR/logs/star_sim.log"

samtools index -@ $THREADS "$OUTDIR/bam/sim_Aligned.sortedByCoord.out.bam"
echo "Step 2 complete."

echo ""
echo "Step 3: Computing CDS uniqueness factors..."

python3 - << PYEOF
import pysam, csv
from collections import defaultdict
from pathlib import Path

outdir   = Path("$OUTDIR")
bam_path = outdir / "bam/sim_Aligned.sortedByCoord.out.bam"
r1_path  = outdir / "simreads/sim_R1.fastq"

sim_counts = defaultdict(int)
with open(r1_path) as fh:
    for line in fh:
        if line.startswith('@'):
            sim_counts[line[1:].split('|')[0]] += 1

unique_back = defaultdict(int)
total_back  = defaultdict(int)
bam = pysam.AlignmentFile(str(bam_path), 'rb')
for read in bam.fetch():
    if read.is_supplementary or read.is_secondary or read.is_unmapped: continue
    cid = read.query_name.split('|')[0]
    total_back[cid] += 1
    if read.mapping_quality == 255: unique_back[cid] += 1
bam.close()

records = []
for cid, n_sim in sim_counts.items():
    n_u = unique_back.get(cid, 0); n_t = total_back.get(cid, 0)
    records.append({
        'cds_id':            cid,
        'n_positions':       n_sim,
        'n_unique_back':     n_u,
        'n_multi_back':      n_t - n_u,
        'n_unmapped':        n_sim - n_t,
        'uniqueness_factor': round(n_u / n_sim, 6) if n_sim > 0 else 0.0,
    })

records.sort(key=lambda r: r['uniqueness_factor'])
out = outdir / "cds_uniqueness_factors_genome_SE_L200bp.tsv"
fields = ['cds_id','n_positions','n_unique_back','n_multi_back','n_unmapped','uniqueness_factor']
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(records)

fs = [r['uniqueness_factor'] for r in records]
n  = len(fs)
print(f"  n={n:,}  mean={sum(fs)/n:.4f}  f=1.00: {sum(1 for f in fs if f==1.0):,}  f=0.00: {sum(1 for f in fs if f==0.0):,}")
print(f"  Saved: {out}")
PYEOF

echo "════════════════════════════════════════════════════════"
echo "Done: $OUTDIR/cds_uniqueness_factors_genome_SE_L200bp.tsv"
echo "════════════════════════════════════════════════════════"
