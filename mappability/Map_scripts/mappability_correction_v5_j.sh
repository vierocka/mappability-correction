#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — v5_j
#
# Reference : MANE ref
# Reads     : exon-flank
# Mode      : SE
# Read len  : 75 bp
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$SCRIPT_DIR/../../Ref"
MANE_GFF="$REF/MANE.GRCh38.v1.5.ensembl_genomic.gff"
GENOME_FA="$REF/GCF_000001405.40_GRCh38.p14_genomic.fna"
MANE_RNA="$REF/MANE.GRCh38.v1.5.ensembl_rna.fna"
DEDUP_FA="$REF/GCF_000001405.40_GRCh38.p14_cds_from_genomic.dedup.fna"
GENOME_IDX="$REF/STAR_index_genome"
CDS_IDX="$REF/STAR_index_dedup_cds"
MANE_IDX="$REF/STAR_index_MANE"
OUTDIR="$SCRIPT_DIR/../results/mappability_MANE_exonflank_SE_L75bp"
THREADS=64
READ_LEN=75
FLANK=50

mkdir -p "$OUTDIR/simreads" "$OUTDIR/bam" "$OUTDIR/logs"

echo "════════════════════════════════════════════════════════"
echo "Mappability — v5_j (MANE ref, exon-flank, SE, 75 bp)"
echo "READ_LEN=$READ_LEN  THREADS=$THREADS"
echo "════════════════════════════════════════════════════════"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0: Build chr→NC chromosome name map
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0: Building chr→NC chromosome name map..."
python3 - << PYEOF
import re
from pathlib import Path

genome_fa = Path("$GENOME_FA")
out_map   = Path("$OUTDIR/chr_name_map.tsv")

ncbi_map = {}
print("  Scanning FASTA headers...", flush=True)
with open(genome_fa) as fh:
    for line in fh:
        if not line.startswith('>'): continue
        name = line[1:].split()[0]
        rest = line.strip()
        m = re.search(r'chromosome (\w+)[,\s]', rest)
        if m and name.startswith('NC_'):
            ncbi_map[f'chr{m.group(1)}'] = name
        if 'mitochondrion' in rest.lower() or 'mitochondrial' in rest.lower():
            ncbi_map['chrM'] = name; ncbi_map['chrMT'] = name

with open(out_map, 'w') as fh:
    for k, v in sorted(ncbi_map.items()):
        fh.write(f"{k}\t{v}\n")
print(f"  Chromosomes mapped: {len(ncbi_map)}")
PYEOF
echo "Step 0 complete."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0b: Build STAR index (MANE RNA) — shared; skipped if already exists
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0b: Checking STAR index (MANE RNA)..."
mkdir -p "$MANE_IDX"
if [ ! -f "$MANE_IDX"/SA ]; then
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$MANE_IDX" \
        --genomeFastaFiles    "$MANE_RNA" \
        --genomeSAindexNbases 12 \
        --genomeChrBinNbits   11 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR/logs/star_index_MANE.log"
    echo "Step 0b complete."
else
    echo "Step 0b: index already exists, skipping."
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Parse GFF, fetch exon ± FLANK from genome, simulate reads
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Simulating exon-flank reads (READ_LEN=$READ_LEN, FLANK=$FLANK, SE)..."

python3 - << PYEOF
import subprocess, re
from pathlib import Path
from collections import defaultdict

gff_path  = Path("$MANE_GFF")
genome_fa = Path("$GENOME_FA")
outdir    = Path("$OUTDIR")
read_len  = int("$READ_LEN")
flank     = int("$FLANK")


chr_map = {}
map_file = outdir / "chr_name_map.tsv"
with open(map_file) as fh:
    for line in fh:
        parts = line.strip().split('\t')
        if len(parts) == 2: chr_map[parts[0]] = parts[1]
print(f"  Loaded {len(chr_map)} chr->NC mappings", flush=True)

tx_exons = defaultdict(list)
with open(gff_path) as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        chrom = chr_map.get(f[0], f[0])
        start = int(f[3]) - 1; end = int(f[4]); strand = f[6]
        tid = None
        for part in f[8].split(';'):
            if part.startswith('transcript_id='):
                tid = part.split('=', 1)[1].strip(); break
        if tid: tx_exons[tid].append((chrom, start, end, strand))
print(f"  Transcripts: {len(tx_exons):,}", flush=True)

fai = Path(str(genome_fa) + '.fai')
if not fai.exists():
    import subprocess as sp
    sp.run(['samtools', 'faidx', str(genome_fa)], check=True)

r1_path = outdir / 'simreads/sim_R1.fastq'

qual = 'I' * read_len
n_reads = n_skip = 0

with open(r1_path, 'w') as r1:
    for tid, exons in tx_exons.items():
        for ex_idx, (chrom, ex_start, ex_end, strand) in enumerate(exons):
            win_start = max(0, ex_start - flank)
            win_end   = ex_end + flank
            region    = f"{chrom}:{win_start + 1}-{win_end}"
            try:
                res = subprocess.run(['samtools', 'faidx', str(genome_fa), region],
                                     capture_output=True, text=True, check=True)
                seq = ''.join(res.stdout.split('\n')[1:]).upper()
            except subprocess.CalledProcessError:
                n_skip += 1; continue
            L = len(seq)
            if L < read_len: continue
            for pos in range(L - read_len + 1):
                fwd = seq[pos : pos + read_len]
                if fwd.count('N') > read_len // 5: continue
                name = f"{tid}|{ex_idx}|{pos}"
                r1.write(f'@{name}\n{fwd}\n+\n{qual}\n')
                n_reads += 1

print(f"  Exon windows skipped: {n_skip}")
print(f"  Total reads simulated: {n_reads:,}", flush=True)
PYEOF
echo "Step 1 complete."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: STAR alignment — SE, MANE ref
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 2: STAR alignment (SE, MANE ref, $THREADS threads)..."

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
import pysam
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
print(f"  Transcripts: {len(sim_counts):,}  Reads: {sum(sim_counts.values()):,}", flush=True)

unique_back = defaultdict(int)
total_back  = defaultdict(int)
bam = pysam.AlignmentFile(str(bam_path), 'rb')
for read in bam.fetch():
    if read.is_supplementary or read.is_secondary or read.is_unmapped: continue
    tid = read.query_name.split('|')[0]
    total_back[tid] += 1
    if read.mapping_quality == 255: unique_back[tid] += 1
bam.close()

import csv
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
out = outdir / "transcript_uniqueness_factors_MANE_exonflank_SE_L75bp.tsv"
fields = ['transcript_id','n_positions','n_unique_back','n_multi_back','n_unmapped','uniqueness_factor']
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(records)

fs = [r['uniqueness_factor'] for r in records]
n  = len(fs)
print(f"\n  ── Distribution ─────────────────────────────────────")
print(f"  n         : {n:,}")
print(f"  mean      : {sum(fs)/n:.4f}")
print(f"  f=1.00    : {sum(1 for f in fs if f==1.0):,}")
print(f"  f>=0.90   : {sum(1 for f in fs if f>=0.90):,}")
print(f"  f 0.10-0.90: {sum(1 for f in fs if 0.10<=f<1.0):,}")
print(f"  f<0.10    : {sum(1 for f in fs if f<0.10):,}")
print(f"  f=0.00    : {sum(1 for f in fs if f==0.0):,}")
print(f"\n  Saved: {out}")
PYEOF

echo ""
echo "════════════════════════════════════════════════════════"
echo "Done: $OUTDIR/transcript_uniqueness_factors_MANE_exonflank_SE_L75bp.tsv"
echo "════════════════════════════════════════════════════════"
