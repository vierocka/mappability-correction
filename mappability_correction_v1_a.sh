#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — genomic exon-flank approach
#
# GFF format confirmed (MANE GRCh38.v1.5.ensembl_genomic.gff):
#   - feature type field [2]: "exon"
#   - attributes field [8]: key=value; no quotes; semicolon-separated
#   - transcript_id=ENST00000376030.7  (with version, matches fC TPM files)
#   - chromosome: "chr1" format
#
# Chromosome naming:
#   GFF uses "chr1" names; the NCBI FASTA uses "NC_000001.11" names.
#   Step 0 builds the chr→NC map; samtools faidx is called with the NC name.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$SCRIPT_DIR/Ref"
MANE_GFF="$REF/MANE.GRCh38.v1.5.ensembl_genomic.gff"
GENOME_FA="$REF/GCF_000001405.40_GRCh38.p14_genomic.fna"
GENOME_IDX="$REF/STAR_index_genome"
OUTDIR="$SCRIPT_DIR/results/mappability_genomic"
THREADS=64
READ_LEN=100   # match actual R1 read length (check with: zcat R1.fastq.gz | sed -n '2p' | wc -c)
FLANK=50       # bp of intronic flanking sequence each side of each exon

# GFF "chr1" names are mapped to NCBI "NC_000001.11" names via Step 0.
CHR_NAMES="ncbi"

mkdir -p "$OUTDIR/simreads" "$OUTDIR/bam" "$OUTDIR/logs"

echo "════════════════════════════════════════════════════════"
echo "MANE mappability — genomic exon-flank approach"
echo "READ_LEN=$READ_LEN  FLANK=$FLANK  THREADS=$THREADS"
echo "GFF chromosome naming: $CHR_NAMES"
echo "════════════════════════════════════════════════════════"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0: Build chr → FASTA sequence-name map (only if CHR_NAMES=ncbi)
# ─────────────────────────────────────────────────────────────────────────────
if [ "$CHR_NAMES" = "ncbi" ]; then
    echo ""
    echo "Step 0: Building chr→NC chromosome name map..."
    python3 - << PYEOF
# Build map from NCBI genome FASTA headers
# Header format: >NC_000001.11 Homo sapiens chromosome 1, GRCh38.p14 primary ...
import gzip, re
from pathlib import Path

genome_fa = Path("$GENOME_FA")
out_map   = Path("$OUTDIR/chr_name_map.tsv")

# NCBI standard: chromosome N → NC_00000N.XX
# chrM → NC_012920.1 (mitochondrial)
# chrX → NC_000023.11
# chrY → NC_000024.10
ncbi_map = {}
opener = gzip.open if str(genome_fa).endswith('.gz') else open
print("  Scanning FASTA headers...", flush=True)
with opener(genome_fa, 'rt') as fh:
    for line in fh:
        if not line.startswith('>'): continue
        line = line.strip()
        name = line[1:].split()[0]   # e.g. NC_000001.11
        rest = line[1:]
        # Look for "chromosome N," or "chromosome X," or "chromosome Y,"
        m = re.search(r'chromosome (\w+)[,\s]', rest)
        if m:
            chrom_num = m.group(1)
            chr_name  = f"chr{chrom_num}"
            ncbi_map[chr_name] = name
        # Mitochondrial
        if 'mitochondrion' in rest.lower() or 'mitochondrial' in rest.lower():
            ncbi_map['chrM'] = name
            ncbi_map['chrMT'] = name

with open(out_map, 'w') as fh:
    for k, v in sorted(ncbi_map.items()):
        fh.write(f"{k}\t{v}\n")

print(f"  Chromosomes mapped: {len(ncbi_map)}")
for k in sorted(ncbi_map)[:5]:
    print(f"    {k} → {ncbi_map[k]}")
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
# STEP 1: Parse GFF, fetch exon ± FLANK from genome, simulate reads
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Parsing GFF and simulating reads..."

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

# ── Load chr→NC map if needed ─────────────────────────────────────────────────
chr_map = {}
if chr_names == "ncbi":
    map_file = outdir / "chr_name_map.tsv"
    with open(map_file) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                chr_map[parts[0]] = parts[1]
    print(f"  Loaded {len(chr_map)} chr→NC mappings", flush=True)

def resolve_chrom(chr_gff):
    return chr_map.get(chr_gff, chr_gff)

# ── Parse GFF for exon features ───────────────────────────────────────────────
print("  Parsing GFF exon features...", flush=True)
tx_exons = defaultdict(list)    # {transcript_id: [(chrom_fasta, start, end, strand)]}

opener = gzip.open if str(gff_path).endswith('.gz') else open
n_exons = 0
with opener(gff_path, 'rt') as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue

        chrom_gff = f[0]
        start     = int(f[3]) - 1   # GFF 1-based → 0-based
        end       = int(f[4])        # GFF end inclusive → Python exclusive
        strand    = f[6]

        # Parse transcript_id from attributes
        # Format: transcript_id=ENST00000376030.7  (no quotes, GFF3)
        tid = None
        for part in f[8].split(';'):
            if part.startswith('transcript_id='):
                tid = part.split('=', 1)[1].strip()
                break

        if tid is None: continue

        chrom_fa = resolve_chrom(chrom_gff)
        tx_exons[tid].append((chrom_fa, start, end, strand))
        n_exons += 1

print(f"  Exon features parsed: {n_exons:,}")
print(f"  Unique transcripts: {len(tx_exons):,}", flush=True)

# ── Index genome FASTA ────────────────────────────────────────────────────────
fai = Path(str(genome_fa) + '.fai')
if not fai.exists():
    print("  Indexing genome with samtools faidx...", flush=True)
    subprocess.run(['samtools', 'faidx', str(genome_fa)], check=True)

# ── Simulate reads from exon ± flank windows ─────────────────────────────────
print("  Simulating reads (step=1, READ_LEN={read_len}, FLANK={flank})...".format(
    read_len=read_len, flank=flank), flush=True)

r1_path = outdir / "simreads/sim_R1.fastq"
r2_path = outdir / "simreads/sim_R2.fastq"
qual    = 'I' * read_len

n_reads = 0
n_skip  = 0

with open(r1_path, 'w') as r1, open(r2_path, 'w') as r2:
    for tid, exons in tx_exons.items():
        for ex_idx, (chrom_fa, ex_start, ex_end, strand) in enumerate(exons):
            win_start = max(0, ex_start - flank)
            win_end   = ex_end + flank
            # samtools faidx: 1-based inclusive
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
                r2.write(f"@{name}\n{revcomp(fwd)}\n+\n{qual}\n")
                n_reads += 1

print(f"  Exon windows skipped (sequence fetch error): {n_skip}")
print(f"  Total reads simulated: {n_reads:,}", flush=True)
PYEOF

echo "Step 1 complete."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: STAR alignment — splice-aware (introns enabled)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 2: STAR alignment (splice-aware, $THREADS threads)..."

STAR \
    --runThreadN      $THREADS \
    --genomeDir       "$GENOME_IDX" \
    --readFilesIn     "$OUTDIR/simreads/sim_R1.fastq" \
                      "$OUTDIR/simreads/sim_R2.fastq" \
    --outSAMtype      BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM \
    --outSAMmultNmax  1 \
    --outFilterMultimapNmax 40 \
    --outBAMsortingThreadN $THREADS \
    --outBAMsortingBinsN   20 \
    --limitBAMsortRAM      160000000000 \
    --outFileNamePrefix "$OUTDIR/bam/sim_" \
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
            # name: tid|ex_idx|pos — split on first two '|', tid may contain '|'
            tid = line[1:].split('|')[0]
            sim_counts[tid] += 1

print(f"  Transcripts: {len(sim_counts):,}  "
      f"Reads: {sum(sim_counts.values()):,}", flush=True)

print("  Counting unique recoveries from BAM...", flush=True)
unique_back = defaultdict(int)
total_back  = defaultdict(int)

bam = pysam.AlignmentFile(str(bam_path), 'rb')
for read in bam.fetch():
    if read.is_supplementary or read.is_secondary or read.is_read2:
        continue
    if read.is_unmapped: continue
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
out = outdir / "transcript_uniqueness_factors_genomic_L100bp.tsv"
df.to_csv(out, sep='\t', index=False)

print(f"\n  ── Distribution ─────────────────────────────────────")
print(df['uniqueness_factor'].describe().round(4).to_string())
print(f"\n  f=1.00 (perfect):   {(df.uniqueness_factor==1.0).sum():,}")
print(f"  f>=0.90:            {(df.uniqueness_factor>=0.90).sum():,}")
print(f"  f 0.10–0.90:        "
      f"{((df.uniqueness_factor>=0.10)&(df.uniqueness_factor<0.90)).sum():,}")
print(f"  f<0.10 (exclude):   {(df.uniqueness_factor<0.10).sum():,}")
print(f"\n  ── Worst 15 ─────────────────────────────────────────")
print(df.head(15)[['transcript_id','n_positions',
                   'n_unique_back','uniqueness_factor']].to_string(index=False))
print(f"\n  Saved: {out}")
PYEOF

echo ""
echo "════════════════════════════════════════════════════════"
echo "Done: $OUTDIR/transcript_uniqueness_factors_genomic_L100bp.tsv"
echo "Run: python3 apply_uniqueness_correction.py"
echo "════════════════════════════════════════════════════════"
