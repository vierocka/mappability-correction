#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — v7_g
#
# Read source : deduplicated CDS (93,088 sequences)
# Reference   : deduplicated CDS (self-mapping)
# Mode        : SE
# Read len    : 100 bp
#
# Reads simulated from dedup CDS and aligned back to the dedup CDS reference.
# Measures within-CDS-set uniqueness: which sequences share enough k-mers that
# even within the deduplicated reference, reads cannot be placed unambiguously.
# Reveals residual redundancy that survived seqkit rmdup (sequences with
# <100% identity but overlapping k-mer content at the chosen read length).
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$SCRIPT_DIR/../../Ref"
DEDUP_FA="$REF/GCF_000001405.40_GRCh38.p14_cds_from_genomic.dedup.fna"
CDS_IDX="$REF/STAR_index_dedup_cds"
OUTDIR="$SCRIPT_DIR/../results/mappability_dedup_cds_reads_selfmap_SE"
THREADS=64
READ_LEN=100

mkdir -p "$OUTDIR/simreads" "$OUTDIR/bam" "$OUTDIR/logs"

echo "════════════════════════════════════════════════════════"
echo "Mappability — v7_g (dedup CDS reads, dedup CDS self-map, SE, 100 bp)"
echo "READ_LEN=$READ_LEN  THREADS=$THREADS"
echo "════════════════════════════════════════════════════════"

echo ""
echo "Step 0: Checking STAR index (dedup CDS)..."
mkdir -p "$CDS_IDX"
if [ ! -f "$CDS_IDX/SA" ]; then
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$CDS_IDX" \
        --genomeFastaFiles    "$DEDUP_FA" \
        --genomeSAindexNbases 12 \
        --genomeChrBinNbits   11 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR/logs/star_index_cds.log"
    echo "Step 0 complete."
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
echo "Step 2: STAR alignment (SE, dedup CDS self-map, $THREADS threads)..."

STAR \
    --runThreadN            $THREADS \
    --genomeDir             "$CDS_IDX" \
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
echo "Step 3: Computing CDS self-mapping uniqueness factors..."

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
print(f"  CDS seqs: {len(sim_counts):,}  Reads: {sum(sim_counts.values()):,}", flush=True)

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
out = outdir / "cds_uniqueness_factors_selfmap_SE.tsv"
fields = ['cds_id','n_positions','n_unique_back','n_multi_back','n_unmapped','uniqueness_factor']
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(records)

fs = [r['uniqueness_factor'] for r in records]
n  = len(fs)
print(f"\n  ── Distribution ─────────────────────────────────────")
print(f"  n          : {n:,}")
print(f"  mean       : {sum(fs)/n:.4f}")
print(f"  f=1.00     : {sum(1 for f in fs if f==1.0):,}")
print(f"  f>=0.90    : {sum(1 for f in fs if f>=0.90):,}")
print(f"  f<0.10     : {sum(1 for f in fs if f<0.10):,}")
print(f"  f=0.00     : {sum(1 for f in fs if f==0.0):,}")
print(f"\n  Saved: {out}")
PYEOF

echo ""
echo "════════════════════════════════════════════════════════"
echo "Done: $OUTDIR/cds_uniqueness_factors_selfmap_SE.tsv"
echo "════════════════════════════════════════════════════════"
