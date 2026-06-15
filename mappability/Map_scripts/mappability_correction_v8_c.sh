#!/bin/bash
# v8_c — MANE RNA reads | proper PE simulation | 150 bp | MANE RNA reference (upper bound)
#
# Reuses reads from v8_a (SIM_DIR). Aligns against MANE RNA index (no pseudogenes).
# UF → 1.0 expected for most transcripts — self-mapping upper bound.
# With proper PE geometry: confirms that the concordance constraint works correctly
# even when insert size matches the reference (no spurious multimapping).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$SCRIPT_DIR/../../Ref"
MANE_RNA="$REF/MANE.GRCh38.v1.5.ensembl_rna.fna"
MANE_IDX="$REF/STAR_index_MANE"
SIM_DIR="$SCRIPT_DIR/../results/v8_simreads_MANE_RNA_PE_L150bp_proper"
OUTDIR="$SCRIPT_DIR/../results/mappability_MANE_RNA_PE_L150bp_proper"
THREADS=64

mkdir -p "$OUTDIR/bam" "$OUTDIR/logs"

echo "════════════════════════════════════════════════════════"
echo "Mappability — v8_c (MANE RNA reads, proper PE, 150 bp, MANE RNA ref)"
echo "Reads from: $SIM_DIR"
echo "════════════════════════════════════════════════════════"

# ── Step 0: STAR index (MANE RNA) ─────────────────────────────────────────────

echo ""
echo "Step 0: Checking STAR MANE RNA index..."
mkdir -p "$MANE_IDX"
if [ ! -f "$MANE_IDX/SA" ]; then
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$MANE_IDX" \
        --genomeFastaFiles    "$MANE_RNA" \
        --genomeSAindexNbases 12 \
        --genomeChrBinNbits   11 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR/logs/star_index.log"
else
    echo "  Index exists, skipping."
fi

# ── Step 1: Check reads ───────────────────────────────────────────────────────

echo ""
echo "Step 1: Checking for reads from v8_a..."
if [ ! -f "$SIM_DIR/sim_R1.fastq" ]; then
    echo "ERROR: reads not found at $SIM_DIR — run v8_a first." >&2
    exit 1
fi
echo "  Reads found, proceeding."

# ── Step 2: STAR alignment ────────────────────────────────────────────────────

echo ""
echo "Step 2: STAR alignment (PE, MANE RNA ref, $THREADS threads)..."

STAR \
    --runThreadN            $THREADS \
    --genomeDir             "$MANE_IDX" \
    --readFilesIn           "$SIM_DIR/sim_R1.fastq" \
                            "$SIM_DIR/sim_R2.fastq" \
    --outSAMtype            BAM SortedByCoordinate \
    --outSAMattributes      NH HI AS NM \
    --outSAMmultNmax        1 \
    --outFilterMultimapNmax 40 \
    --alignIntronMax        1 \
    --outBAMsortingThreadN  $THREADS \
    --outBAMsortingBinsN    20 \
    --limitBAMsortRAM       160000000000 \
    --outFileNamePrefix     "$OUTDIR/bam/sim_" \
    2>&1 | tee "$OUTDIR/logs/star_sim.log"

samtools index -@ $THREADS "$OUTDIR/bam/sim_Aligned.sortedByCoord.out.bam"
echo "Step 2 complete."

# ── Step 3: Uniqueness factors ────────────────────────────────────────────────

echo ""
echo "Step 3: Computing transcript uniqueness factors..."

python3 - << PYEOF
import pysam, csv
from collections import defaultdict
from pathlib import Path

simdir   = Path("$SIM_DIR")
outdir   = Path("$OUTDIR")
bam_path = outdir / "bam/sim_Aligned.sortedByCoord.out.bam"
r1_path  = simdir / "sim_R1.fastq"

sim_counts = defaultdict(int)
with open(r1_path) as fh:
    for line in fh:
        if line.startswith('@'):
            sim_counts[line[1:].split('|')[0]] += 1

unique_back = defaultdict(int)
total_back  = defaultdict(int)
bam = pysam.AlignmentFile(str(bam_path), 'rb')
for read in bam.fetch():
    if read.is_supplementary or read.is_secondary or read.is_read2: continue
    tid = read.query_name.split('|')[0]
    total_back[tid] += 1
    if read.mapping_quality == 255: unique_back[tid] += 1
bam.close()

records = []
for tid, n_sim in sim_counts.items():
    n_u = unique_back.get(tid, 0); n_t = total_back.get(tid, 0)
    records.append({
        'transcript_id':     tid,
        'n_positions':       n_sim,
        'n_unique_back':     n_u,
        'n_multi_back':      n_t - n_u,
        'n_unmapped':        n_sim - n_t,
        'uniqueness_factor': round(n_u / n_sim, 6) if n_sim > 0 else 0.0,
    })

records.sort(key=lambda r: r['uniqueness_factor'])
out = outdir / "transcript_uniqueness_factors_MANE_RNA_PE_L150bp_proper.tsv"
fields = ['transcript_id','n_positions','n_unique_back','n_multi_back','n_unmapped','uniqueness_factor']
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(records)

fs = [r['uniqueness_factor'] for r in records]
n  = len(fs)
print(f"  n={n:,}  mean={sum(fs)/n:.4f}  UF=1.0: {sum(1 for f in fs if f==1.0):,}  UF=0.0: {sum(1 for f in fs if f==0.0):,}")
print(f"  Saved: {out}")
PYEOF

echo "════════════════════════════════════════════════════════"
echo "Done: $OUTDIR/transcript_uniqueness_factors_MANE_RNA_PE_L150bp_proper.tsv"
echo "════════════════════════════════════════════════════════"
