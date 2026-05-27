#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — v1_b
#
# Identical to v1 (genomic exon-flank approach, full genome reference)
# EXCEPT: single-end alignment (R2 removed).
#
# Comparison axis v1 → v1_b: isolates the effect of PE vs SE on full genome.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$SCRIPT_DIR/../Ref"
MANE_GFF="$REF/MANE.GRCh38.v1.5.ensembl_genomic.gff"
GENOME_FA="$REF/GCF_000001405.40_GRCh38.p14_genomic.fna"
GENOME_IDX="$REF/STAR_index_genome"
OUTDIR="$SCRIPT_DIR/results/mappability_genomic_SE"
THREADS=64
READ_LEN=100
FLANK=50
CHR_NAMES="ncbi"

mkdir -p "$OUTDIR/simreads" "$OUTDIR/bam" "$OUTDIR/logs"

echo "════════════════════════════════════════════════════════"
echo "MANE mappability — v1_b (full genome, exon-flank, SE)"
echo "READ_LEN=$READ_LEN  FLANK=$FLANK  THREADS=$THREADS"
echo "════════════════════════════════════════════════════════"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0: Build chr → FASTA sequence-name map (identical to v1)
# ─────────────────────────────────────────────────────────────────────────────
if [ "$CHR_NAMES" = "ncbi" ]; then
    echo ""
    echo "Step 0: Building chr→NC chromosome name map..."
    python3 - << PYEOF
import gzip, re
from pathlib import Path

genome_fa = Path("$GENOME_FA")
out_map   = Path("$OUTDIR/chr_name_map.tsv")

ncbi_map = {}
opener = gzip.open if str(genome_fa).endswith('.gz') else open
print("  Scanning FASTA headers...", flush=True)
with opener(genome_fa, 'rt') as fh:
    for line in fh:
        if not line.startswith('>'): continue
        line = line.strip()
        name = line[1:].split()[0]
        rest = line[1:]
        m = re.search(r'chromosome (\w+)[,\s]', rest)
        if m:
            ncbi_map[f"chr{m.group(1)}"] = name
        if 'mitochondrion' in rest.lower() or 'mitochondrial' in rest.lower():
            ncbi_map['chrM']  = name
            ncbi_map['chrMT'] = name

with open(out_map, 'w') as fh:
    for k, v in sorted(ncbi_map.items()):
        fh.write(f"{k}\t{v}\n")

print(f"  Chromosomes mapped: {len(ncbi_map)}")
PYEOF
    echo "Step 0 complete."
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0b: Build STAR index from full genome (once; shared by all genome scripts)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0b: Building STAR index from full genome..."
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
    echo "Step 0b complete."
else
    echo "Step 0b: Index already exists, skipping."
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Parse GFF, fetch exon ± FLANK from genome, simulate reads (SE only)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Parsing GFF and simulating reads (single-end)..."

python3 - << PYEOF
import gzip, subprocess
from pathlib import Path
from collections import defaultdict

gff_path  = Path("$MANE_GFF")
genome_fa = Path("$GENOME_FA")
outdir    = Path("$OUTDIR")
read_len  = $READ_LEN
flank     = $FLANK
chr_names = "$CHR_NAMES"

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def revcomp(s): return s.translate(RC)[::-1]

chr_map = {}
if chr_names == "ncbi":
    with open(outdir / "chr_name_map.tsv") as fh:
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                chr_map[parts[0]] = parts[1]
    print(f"  Loaded {len(chr_map)} chr→NC mappings", flush=True)

def resolve_chrom(c):
    return chr_map.get(c, c)

print("  Parsing GFF exon features...", flush=True)
tx_exons = defaultdict(list)
opener = gzip.open if str(gff_path).endswith('.gz') else open
n_exons = 0
with opener(gff_path, 'rt') as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        tid = None
        for part in f[8].split(';'):
            if part.startswith('transcript_id='):
                tid = part.split('=', 1)[1].strip()
                break
        if tid is None: continue
        tx_exons[tid].append((resolve_chrom(f[0]), int(f[3])-1, int(f[4]), f[6]))
        n_exons += 1

print(f"  Exon features parsed: {n_exons:,}")
print(f"  Unique transcripts: {len(tx_exons):,}", flush=True)

fai = Path(str(genome_fa) + '.fai')
if not fai.exists():
    print("  Indexing genome with samtools faidx...", flush=True)
    subprocess.run(['samtools', 'faidx', str(genome_fa)], check=True)

print(f"  Simulating reads (SE, step=1, READ_LEN={read_len}, FLANK={flank})...", flush=True)

r1_path = outdir / "simreads/sim_R1.fastq"
qual    = 'I' * read_len
n_reads = 0
n_skip  = 0

with open(r1_path, 'w') as r1:
    for tid, exons in tx_exons.items():
        for ex_idx, (chrom_fa, ex_start, ex_end, strand) in enumerate(exons):
            win_start = max(0, ex_start - flank)
            win_end   = ex_end + flank
            region = f"{chrom_fa}:{win_start + 1}-{win_end}"
            try:
                res = subprocess.run(
                    ['samtools', 'faidx', str(genome_fa), region],
                    capture_output=True, text=True, check=True)
                seq = ''.join(res.stdout.split('\n')[1:]).upper()
            except subprocess.CalledProcessError:
                n_skip += 1
                continue
            L = len(seq)
            if L < read_len:
                continue
            for pos in range(L - read_len + 1):
                fwd = seq[pos : pos + read_len]
                if fwd.count('N') > read_len // 5:
                    continue
                name = f"{tid}|{ex_idx}|{pos}"
                r1.write(f"@{name}\n{fwd}\n+\n{qual}\n")
                n_reads += 1

print(f"  Exon windows skipped: {n_skip}")
print(f"  Total reads simulated: {n_reads:,}", flush=True)
PYEOF

echo "Step 1 complete."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: STAR alignment — single-end, splice-aware
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 2: STAR alignment single-end ($THREADS threads)..."

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
out = outdir / "transcript_uniqueness_factors_genomic_SE_L100bp.tsv"
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
echo "Done: $OUTDIR/transcript_uniqueness_factors_genomic_SE_L100bp.tsv"
echo "════════════════════════════════════════════════════════"
