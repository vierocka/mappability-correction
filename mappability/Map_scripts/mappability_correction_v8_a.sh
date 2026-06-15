#!/bin/bash
# v8_a — MANE RNA reads | proper PE simulation | 150 bp | full genome
#
# Simulates realistic PE 150 bp reads: fragment = 300 bp, R1 from the 5' end,
# R2 (reverse complement) from the 3' end. Unlike v4_d, R1 and R2 come from
# opposite ends of the same fragment — matching real Illumina PE library geometry.
# Compare to v4_d (same parameters, fake PE: R2 = revcomp(R1) at the same window).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$SCRIPT_DIR/../../Ref"
MANE_RNA="$REF/MANE.GRCh38.v1.5.ensembl_rna.fna"
GENOME_FA="$REF/GCF_000001405.40_GRCh38.p14_genomic.fna"
GENOME_IDX="$REF/STAR_index_genome"
SIM_DIR="$SCRIPT_DIR/../results/v8_simreads_MANE_RNA_PE_L150bp_proper"
OUTDIR="$SCRIPT_DIR/../results/mappability_genomic_RNA_PE_L150bp_proper"
THREADS=64
READ_LEN=150

mkdir -p "$SIM_DIR" "$OUTDIR/bam" "$OUTDIR/logs"

echo "════════════════════════════════════════════════════════"
echo "Mappability — v8_a (MANE RNA reads, proper PE, 150 bp, full genome)"
echo "READ_LEN=$READ_LEN  FRAG_LEN=$((READ_LEN * 2))  THREADS=$THREADS"
echo "════════════════════════════════════════════════════════"

# ── Step 0: STAR index ────────────────────────────────────────────────────────

echo ""
echo "Step 0: Checking STAR genome index..."
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
    echo "  Index exists, skipping."
fi

# ── Step 1: Simulate reads ────────────────────────────────────────────────────

echo ""
echo "Step 1: Simulating MANE RNA reads (proper PE, 150 bp, 300 bp fragment)..."

if [ -f "$SIM_DIR/sim_R1.fastq" ]; then
    echo "  Reads exist ($SIM_DIR), skipping."
else

python3 - << PYEOF
from pathlib import Path

mane_rna = Path("$MANE_RNA")
simdir   = Path("$SIM_DIR")
read_len = int("$READ_LEN")
frag_len = read_len * 2  # 300 bp: R1 tiles first half, R2 tiles second half

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
            rev = seq[pos + read_len : pos + frag_len]   # second half of fragment
            if fwd.count('N') > read_len // 5: continue
            if rev.count('N') > read_len // 5: continue
            r1.write(f'@{seq_id}|{pos}\n{fwd}\n+\n{qual}\n')
            r2.write(f'@{seq_id}|{pos}\n{revcomp(rev)}\n+\n{qual}\n')
            n_reads += 1

print(f"  Transcripts: {n_seq:,}  (too short for {frag_len} bp fragment: {n_skip})")
print(f"  Read pairs simulated: {n_reads:,}", flush=True)
PYEOF

fi
echo "Step 1 complete."

# ── Step 2: STAR alignment ────────────────────────────────────────────────────

echo ""
echo "Step 2: STAR alignment (PE, full genome, $THREADS threads)..."

STAR \
    --runThreadN            $THREADS \
    --genomeDir             "$GENOME_IDX" \
    --readFilesIn           "$SIM_DIR/sim_R1.fastq" \
                            "$SIM_DIR/sim_R2.fastq" \
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
out = outdir / "transcript_uniqueness_factors_genomic_RNA_PE_L150bp_proper.tsv"
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
echo "Done: $OUTDIR/transcript_uniqueness_factors_genomic_RNA_PE_L150bp_proper.tsv"
echo "════════════════════════════════════════════════════════"
