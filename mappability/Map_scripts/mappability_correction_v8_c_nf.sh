#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — v8_c_nf — NF REFERENCE
#
# Reuses reads simulated by mappability_correction_v8_a_nf.sh (MANE RNA,
# proper PE, 150 bp) — run v8_a_nf first. Aligns against MANE RNA (nf, self-map)
# instead of the nf genome, to isolate the competitor-set effect under
# realistic PE geometry (same comparison as the old v8_b/v8_c).
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS="$SCRIPT_DIR/../results"

REF="$SCRIPT_DIR/../../Ref"
GENOME_FA_NF="$REF/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
GTF_NF="$REF/Homo_sapiens.GRCh38.116.gtf"
MANE_RNA_NF="$REF/MANE_RNA_nf.fna"
MANE_IDX_NF="$REF/STAR_index_MANE_nf"
DEDUP_CDS_NF_RAW="$REF/Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_raw.fna"
DEDUP_CDS_NF="$REF/Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna"
CDS_IDX_NF="$REF/STAR_index_dedup_cds_nf"

SIM_DIR="$RESULTS/mappability_genomic_RNA_PE_L150bp_proper_nf/simreads"
OUTDIR="$RESULTS/mappability_MANE_RNA_PE_L150bp_proper_nf"
THREADS=64

mkdir -p "$OUTDIR/bam" "$OUTDIR/logs"

module purge
module load bio/STAR/2.7.11b-GCC-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0
module load lang/Python/3.11.5-GCCcore-13.2.0

ulimit -n 65536

echo "════════════════════════════════════════════════════════"
echo "Mappability — v8_c_nf (MANE RNA reads, proper PE, 150 bp, MANE RNA (nf, self-map))"
echo "Reads from: $SIM_DIR"
echo "════════════════════════════════════════════════════════"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0: Build STAR index (MANE RNA, nf) — shared; skipped if already exists
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0: Checking STAR MANE-RNA-nf index..."
mkdir -p "$MANE_IDX_NF"
if [ ! -f "$MANE_IDX_NF/SA" ]; then
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$MANE_IDX_NF" \
        --genomeFastaFiles    "$MANE_RNA_NF" \
        --genomeSAindexNbases 12 \
        --genomeChrBinNbits   11 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR/logs/star_index_mane_nf.log"
else
    echo "  Index exists, skipping."
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Check reads from v8_a_nf
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Checking for reads from v8_a_nf..."
if [ ! -f "$SIM_DIR/sim_R1.fastq" ]; then
    echo "ERROR: reads not found at $SIM_DIR — run v8_a_nf first." >&2
    exit 1
fi
echo "  Reads found, proceeding."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: STAR alignment (PE, MANE RNA (nf, self-map))
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 2: STAR alignment (PE, MANE RNA (nf, self-map), $THREADS threads)..."

STAR \
    --runThreadN            $THREADS \
    --genomeDir             "$MANE_IDX_NF" \
    --readFilesIn           "$SIM_DIR/sim_R1.fastq" \
                            "$SIM_DIR/sim_R2.fastq" \
    --outSAMtype            BAM SortedByCoordinate \
    --outSAMattributes      NH HI AS NM \
    --outSAMmultNmax        1 \
    --outFilterMultimapNmax 20 \
    --alignIntronMax        1 \
    --outBAMsortingThreadN  $THREADS \
    --outBAMsortingBinsN    20 \
    --limitBAMsortRAM       160000000000 \
    --outFileNamePrefix     "$OUTDIR/bam/sim_" \
    2>&1 | tee "$OUTDIR/logs/star_sim.log"

samtools index -@ $THREADS "$OUTDIR/bam/sim_Aligned.sortedByCoord.out.bam"
echo "Step 2 complete."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Compute transcript uniqueness factors
# ─────────────────────────────────────────────────────────────────────────────
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
out = outdir / "transcript_uniqueness_factors_MANE_RNA_PE_L150bp_proper_nf.tsv"
fields = ['transcript_id','n_positions','n_unique_back','n_multi_back','n_unmapped','uniqueness_factor']
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(records)

fs = [r['uniqueness_factor'] for r in records]; n = len(fs)
print(f"  n={n:,}  mean={sum(fs)/n:.4f}  UF=1.0: {sum(1 for f in fs if f==1.0):,}  UF=0.0: {sum(1 for f in fs if f==0.0):,}")
print(f"  Saved: {out}")
PYEOF

echo "════════════════════════════════════════════════════════"
echo "Done: $OUTDIR/transcript_uniqueness_factors_MANE_RNA_PE_L150bp_proper_nf.tsv"
echo "════════════════════════════════════════════════════════"
