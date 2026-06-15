#!/bin/bash
#SBATCH --cpus-per-task=64
#SBATCH --mem=192gb
#SBATCH --time=48:00:00
#SBATCH --output=mappability_v8_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL

# ─────────────────────────────────────────────────────────────────────────────
# v8 — MANE RNA reads, proper PE 150 bp simulation — HPC submission script
#
# Runs v8_a, v8_b, v8_c sequentially:
#   v8_a: reads vs full genome       → transcript_uniqueness_factors_genomic_RNA_PE_L150bp_proper.tsv
#   v8_b: reads vs dedup CDS         → transcript_uniqueness_factors_dedup_cds_RNA_PE_L150bp_proper.tsv
#   v8_c: reads vs MANE RNA (self)   → transcript_uniqueness_factors_MANE_RNA_PE_L150bp_proper.tsv
#
# PE simulation: fragment = 300 bp; R1 = first 150 bp (forward);
# R2 = last 150 bp (reverse complement). Reads simulated once, reused for all
# three alignments. This reflects real Illumina PE library geometry.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

BASE="./mappability-correction"
REF="$BASE/Ref"
RESULTS="$BASE/mappability/results"
THREADS=64
READ_LEN=150
FRAG_LEN=300

MANE_RNA="$REF/MANE.GRCh38.v1.5.ensembl_rna.fna"
GENOME_FA="$REF/GCF_000001405.40_GRCh38.p14_genomic.fna"
DEDUP_FA="$REF/GCF_000001405.40_GRCh38.p14_cds_from_genomic.dedup.fna"

GENOME_IDX="$REF/STAR_index_genome"
CDS_IDX="$REF/STAR_index_dedup_cds"
MANE_IDX="$REF/STAR_index_MANE"

SIM_DIR="$RESULTS/v8_simreads_MANE_RNA_PE_L150bp_proper"
OUTDIR_A="$RESULTS/mappability_genomic_RNA_PE_L150bp_proper"
OUTDIR_B="$RESULTS/mappability_dedup_cds_RNA_PE_L150bp_proper"
OUTDIR_C="$RESULTS/mappability_MANE_RNA_PE_L150bp_proper"

echo "════════════════════════════════════════════════════════"
echo "v8 mappability — proper PE 150 bp — HPC"
echo "Started : $(date)"
echo "Host    : $(hostname)"
echo "Base    : $BASE"
echo "════════════════════════════════════════════════════════"

module purge
module load bio/STAR/2.7.11b-GCC-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0
module load lang/Python/3.11.5-GCCcore-13.2.0

ulimit -n 65536

mkdir -p "$SIM_DIR" \
         "$OUTDIR_A/bam" "$OUTDIR_A/logs" \
         "$OUTDIR_B/bam" "$OUTDIR_B/logs" \
         "$OUTDIR_C/bam" "$OUTDIR_C/logs"

# ── Step 0: Build STAR indices ────────────────────────────────────────────────

echo ""
echo "Step 0: Building / checking STAR indices..."

mkdir -p "$GENOME_IDX"
if [ ! -f "$GENOME_IDX/SA" ]; then
    echo "  Building full genome index..."
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$GENOME_IDX" \
        --genomeFastaFiles    "$GENOME_FA" \
        --genomeSAindexNbases 14 \
        --genomeChrBinNbits   14 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR_A/logs/star_index_genome.log"
else
    echo "  Full genome index exists."
fi

mkdir -p "$CDS_IDX"
if [ ! -f "$CDS_IDX/SA" ]; then
    echo "  Building dedup CDS index..."
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$CDS_IDX" \
        --genomeFastaFiles    "$DEDUP_FA" \
        --genomeSAindexNbases 12 \
        --genomeChrBinNbits   11 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR_B/logs/star_index_cds.log"
else
    echo "  Dedup CDS index exists."
fi

mkdir -p "$MANE_IDX"
if [ ! -f "$MANE_IDX/SA" ]; then
    echo "  Building MANE RNA index..."
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$MANE_IDX" \
        --genomeFastaFiles    "$MANE_RNA" \
        --genomeSAindexNbases 12 \
        --genomeChrBinNbits   11 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR_C/logs/star_index_mane.log"
else
    echo "  MANE RNA index exists."
fi

# ── Step 1: Simulate reads ────────────────────────────────────────────────────

echo ""
echo "Step 1: Simulating MANE RNA reads (proper PE, ${READ_LEN} bp, ${FRAG_LEN} bp fragment)..."

if [ -f "$SIM_DIR/sim_R1.fastq" ]; then
    echo "  Reads exist, skipping simulation."
else

python3 - << PYEOF
from pathlib import Path

mane_rna = Path("$MANE_RNA")
simdir   = Path("$SIM_DIR")
read_len = int("$READ_LEN")
frag_len = int("$FRAG_LEN")

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def revcomp(s): return s.translate(RC)[::-1]

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

qual    = 'I' * read_len
n_reads = n_skip = n_seq = 0

with open(simdir / 'sim_R1.fastq', 'w') as r1, \
     open(simdir / 'sim_R2.fastq', 'w') as r2:
    for seq_id, seq in parse_fasta(mane_rna):
        L = len(seq)
        if L < frag_len: n_skip += 1; continue
        n_seq += 1
        for pos in range(L - frag_len + 1):
            fwd = seq[pos            : pos + read_len]
            rev = seq[pos + read_len : pos + frag_len]
            if fwd.count('N') > read_len // 5: continue
            if rev.count('N') > read_len // 5: continue
            r1.write(f'@{seq_id}|{pos}\n{fwd}\n+\n{qual}\n')
            r2.write(f'@{seq_id}|{pos}\n{revcomp(rev)}\n+\n{qual}\n')
            n_reads += 1

print(f"  Transcripts: {n_seq:,}  (too short for {frag_len} bp: {n_skip})")
print(f"  Read pairs: {n_reads:,}", flush=True)
PYEOF

fi
echo "Step 1 complete."

# ── Step 2a: STAR → full genome ───────────────────────────────────────────────

echo ""
echo "━━━━ v8_a : MANE RNA reads → full genome ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

STAR \
    --runThreadN            $THREADS \
    --genomeDir             "$GENOME_IDX" \
    --readFilesIn           "$SIM_DIR/sim_R1.fastq" "$SIM_DIR/sim_R2.fastq" \
    --outSAMtype            BAM SortedByCoordinate \
    --outSAMattributes      NH HI AS NM \
    --outSAMmultNmax        1 \
    --outFilterMultimapNmax 40 \
    --outBAMsortingThreadN  $THREADS \
    --outBAMsortingBinsN    20 \
    --limitBAMsortRAM       160000000000 \
    --outFileNamePrefix     "$OUTDIR_A/bam/sim_" \
    2>&1 | tee "$OUTDIR_A/logs/star_sim.log"

samtools index -@ $THREADS "$OUTDIR_A/bam/sim_Aligned.sortedByCoord.out.bam"

python3 - << PYEOF
import pysam, csv
from collections import defaultdict
from pathlib import Path

simdir   = Path("$SIM_DIR")
outdir   = Path("$OUTDIR_A")
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
out = outdir / "transcript_uniqueness_factors_genomic_RNA_PE_L150bp_proper.tsv"
fields = ['transcript_id','n_positions','n_unique_back','n_multi_back','n_unmapped','uniqueness_factor']
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(records)
fs = [r['uniqueness_factor'] for r in records]; n = len(fs)
print(f"  v8_a  n={n:,}  mean={sum(fs)/n:.4f}  UF=1.0: {sum(1 for f in fs if f==1.0):,}  UF=0.0: {sum(1 for f in fs if f==0.0):,}")
print(f"  Saved: {out}", flush=True)
PYEOF

echo "v8_a complete."

# ── Step 2b: STAR → dedup CDS ─────────────────────────────────────────────────

echo ""
echo "━━━━ v8_b : MANE RNA reads → dedup CDS ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

STAR \
    --runThreadN            $THREADS \
    --genomeDir             "$CDS_IDX" \
    --readFilesIn           "$SIM_DIR/sim_R1.fastq" "$SIM_DIR/sim_R2.fastq" \
    --outSAMtype            BAM SortedByCoordinate \
    --outSAMattributes      NH HI AS NM \
    --outSAMmultNmax        1 \
    --outFilterMultimapNmax 40 \
    --alignIntronMax        1 \
    --outBAMsortingThreadN  $THREADS \
    --outBAMsortingBinsN    20 \
    --limitBAMsortRAM       160000000000 \
    --outFileNamePrefix     "$OUTDIR_B/bam/sim_" \
    2>&1 | tee "$OUTDIR_B/logs/star_sim.log"

samtools index -@ $THREADS "$OUTDIR_B/bam/sim_Aligned.sortedByCoord.out.bam"

python3 - << PYEOF
import pysam, csv
from collections import defaultdict
from pathlib import Path

simdir   = Path("$SIM_DIR")
outdir   = Path("$OUTDIR_B")
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
out = outdir / "transcript_uniqueness_factors_dedup_cds_RNA_PE_L150bp_proper.tsv"
fields = ['transcript_id','n_positions','n_unique_back','n_multi_back','n_unmapped','uniqueness_factor']
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(records)
fs = [r['uniqueness_factor'] for r in records]; n = len(fs)
print(f"  v8_b  n={n:,}  mean={sum(fs)/n:.4f}  UF=1.0: {sum(1 for f in fs if f==1.0):,}  UF=0.0: {sum(1 for f in fs if f==0.0):,}")
print(f"  Saved: {out}", flush=True)
PYEOF

echo "v8_b complete."

# ── Step 2c: STAR → MANE RNA (upper bound) ────────────────────────────────────

echo ""
echo "━━━━ v8_c : MANE RNA reads → MANE RNA ref (upper bound) ━━━━━━━━━━━━━━━━━━"

STAR \
    --runThreadN            $THREADS \
    --genomeDir             "$MANE_IDX" \
    --readFilesIn           "$SIM_DIR/sim_R1.fastq" "$SIM_DIR/sim_R2.fastq" \
    --outSAMtype            BAM SortedByCoordinate \
    --outSAMattributes      NH HI AS NM \
    --outSAMmultNmax        1 \
    --outFilterMultimapNmax 40 \
    --alignIntronMax        1 \
    --outBAMsortingThreadN  $THREADS \
    --outBAMsortingBinsN    20 \
    --limitBAMsortRAM       160000000000 \
    --outFileNamePrefix     "$OUTDIR_C/bam/sim_" \
    2>&1 | tee "$OUTDIR_C/logs/star_sim.log"

samtools index -@ $THREADS "$OUTDIR_C/bam/sim_Aligned.sortedByCoord.out.bam"

python3 - << PYEOF
import pysam, csv
from collections import defaultdict
from pathlib import Path

simdir   = Path("$SIM_DIR")
outdir   = Path("$OUTDIR_C")
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
fs = [r['uniqueness_factor'] for r in records]; n = len(fs)
print(f"  v8_c  n={n:,}  mean={sum(fs)/n:.4f}  UF=1.0: {sum(1 for f in fs if f==1.0):,}  UF=0.0: {sum(1 for f in fs if f==0.0):,}")
print(f"  Saved: {out}", flush=True)
PYEOF

echo "v8_c complete."

# ── Summary ───────────────────────────────────────────────────────────────────

echo ""
echo "════════════════════════════════════════════════════════"
echo "All v8 runs complete: $(date)"
echo ""
echo "TSV outputs:"
find "$RESULTS" -name "transcript_uniqueness_factors_*_PE_L150bp_proper.tsv" | sort
echo "════════════════════════════════════════════════════════"
