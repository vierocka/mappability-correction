#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — v5_a
#
# MANE transcriptome reference, single-end alignment.
# Reads simulated from MANE spliced RNA sequences (no intronic flanks).
# Splice junction detection disabled: the reference is already spliced.
#
# Comparison axes:
#   v1_c → v5_a: effect of reference (full genome vs. MANE transcriptome)
#   v5_a → v5_b: effect of read source (MANE RNA vs. genomic exon-flank)
#   v5_a → v5_c: effect of read length (100 vs. 75 bp)
#   v5_a → v5_d: effect of read length (100 vs. 200 bp)
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$SCRIPT_DIR/../../Ref"
MANE_RNA="$REF/MANE.GRCh38.v1.5.ensembl_rna.fna"
MANE_IDX="$REF/STAR_index_MANE"
OUTDIR="$SCRIPT_DIR/../results/mappability_MANE_RNA_SE"
THREADS=64
READ_LEN=100

mkdir -p "$OUTDIR/simreads" "$OUTDIR/bam" "$OUTDIR/logs"

echo "════════════════════════════════════════════════════════"
echo "MANE mappability — v5_a (MANE reference, MANE RNA, SE)"
echo "READ_LEN=$READ_LEN  THREADS=$THREADS"
echo "════════════════════════════════════════════════════════"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0: Build STAR index from MANE RNA FASTA (shared by all v5 scripts)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0: Building STAR index from MANE RNA FASTA..."
mkdir -p "$MANE_IDX"
if [ ! -f "$MANE_IDX/SA" ]; then
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$MANE_IDX" \
        --genomeFastaFiles    "$MANE_RNA" \
        --genomeSAindexNbases 12 \
        --genomeChrBinNbits   11 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR/logs/star_index_MANE.log"
    echo "Step 0 complete."
else
    echo "Step 0: MANE index already exists, skipping."
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Simulate reads from MANE spliced RNA sequences
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Simulating reads from MANE RNA FASTA..."

python3 - << PYEOF
from pathlib import Path

mane_rna = Path("$MANE_RNA")
outdir   = Path("$OUTDIR")
read_len = $READ_LEN

r1_path = outdir / "simreads/sim_R1.fastq"
qual    = 'I' * read_len
n_reads = 0
n_tx    = 0
n_skip  = 0

def parse_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name: yield name, ''.join(seq)
                name = line[1:].split()[0]
                seq  = []
            else:
                seq.append(line.upper())
    if name: yield name, ''.join(seq)

print("  Reading MANE RNA FASTA...", flush=True)
with open(r1_path, 'w') as r1:
    for tid, seq in parse_fasta(mane_rna):
        L = len(seq)
        if L < read_len:
            n_skip += 1
            continue
        n_tx += 1
        for pos in range(L - read_len + 1):
            fwd = seq[pos : pos + read_len]
            if fwd.count('N') > read_len // 5:
                continue
            r1.write(f"@{tid}|{pos}\n{fwd}\n+\n{qual}\n")
            n_reads += 1

print(f"  Transcripts processed: {n_tx:,}  (skipped too short: {n_skip})")
print(f"  Total reads simulated: {n_reads:,}", flush=True)
PYEOF

echo "Step 1 complete."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: STAR alignment — single-end, MANE transcriptome index
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 2: STAR alignment single-end, MANE reference ($THREADS threads)..."

STAR \
    --runThreadN            $THREADS \
    --genomeDir             "$MANE_IDX" \
    --readFilesIn           "$OUTDIR/simreads/sim_R1.fastq" \
    --outSAMtype            BAM SortedByCoordinate \
    --outSAMattributes      NH HI AS NM \
    --outSAMmultNmax        1 \
    --outFilterMultimapNmax 40 \
    --alignIntronMax        1 \
    --alignEndsType         EndToEnd \
    --outBAMsortingThreadN  $THREADS \
    --outBAMsortingBinsN    20 \
    --limitBAMsortRAM       160000000000 \
    --outFileNamePrefix     "$OUTDIR/bam/sim_" \
    2>&1 | tee "$OUTDIR/logs/star_sim.log"

samtools index -@ $THREADS "$OUTDIR/bam/sim_Aligned.sortedByCoord.out.bam"
echo "Step 2 complete."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Compute per-transcript uniqueness factor
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 3: Computing uniqueness factors..."

python3 - << PYEOF
import pysam, pandas as pd
from collections import defaultdict
from pathlib import Path

outdir   = Path("$OUTDIR")
bam_path = outdir / "bam/sim_Aligned.sortedByCoord.out.bam"
r1_path  = outdir / "simreads/sim_R1.fastq"

print("  Counting simulated reads per transcript...", flush=True)
sim_counts = defaultdict(int)
with open(r1_path) as fh:
    for line in fh:
        if line.startswith('@'):
            sim_counts[line[1:].split('|')[0]] += 1

print(f"  Transcripts: {len(sim_counts):,}  Reads: {sum(sim_counts.values()):,}", flush=True)

print("  Counting unique recoveries from BAM...", flush=True)
unique_back = defaultdict(int)
total_back  = defaultdict(int)

bam = pysam.AlignmentFile(str(bam_path), 'rb')
for read in bam.fetch():
    if read.is_supplementary or read.is_secondary or read.is_unmapped:
        continue
    tid = read.query_name.split('|')[0]
    total_back[tid] += 1
    if read.mapping_quality == 255:
        unique_back[tid] += 1
bam.close()

records = []
for tid, n_sim in sim_counts.items():
    n_u = unique_back.get(tid, 0)
    n_t = total_back.get(tid, 0)
    records.append({
        'transcript_id'    : tid,
        'n_positions'      : n_sim,
        'n_unique_back'    : n_u,
        'n_multi_back'     : n_t - n_u,
        'n_unmapped'       : n_sim - n_t,
        'uniqueness_factor': round(n_u / n_sim, 6) if n_sim > 0 else 0.0,
    })

df = pd.DataFrame(records).sort_values('uniqueness_factor')
out = outdir / "transcript_uniqueness_factors_MANE_RNA_SE_L100bp.tsv"
df.to_csv(out, sep='\t', index=False)

print(f"\n  ── Distribution ─────────────────────────────────────")
print(df['uniqueness_factor'].describe().round(4).to_string())
print(f"\n  f=1.00: {(df.uniqueness_factor==1.0).sum():,}")
print(f"  f>=0.90: {(df.uniqueness_factor>=0.90).sum():,}")
print(f"  f<0.10: {(df.uniqueness_factor<0.10).sum():,}")
print(f"  f=0.00: {(df.uniqueness_factor==0.0).sum():,}")
print(f"\n  ── Worst 15 ─────────────────────────────────────────")
print(df.head(15)[['transcript_id','n_positions','n_unique_back','uniqueness_factor']].to_string(index=False))
print(f"\n  Saved: {out}")
PYEOF

echo ""
echo "════════════════════════════════════════════════════════"
echo "Done: $OUTDIR/transcript_uniqueness_factors_MANE_RNA_SE_L100bp.tsv"
echo "════════════════════════════════════════════════════════"
