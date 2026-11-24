#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — v8_b_nf — NF REFERENCE
#
# Reuses reads simulated by mappability_correction_v8_a_nf.sh (MANE RNA,
# proper PE, 150 bp) — run v8_a_nf first. Aligns against dedup CDS (nf)
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
OUTDIR="$RESULTS/mappability_dedup_cds_RNA_PE_L150bp_proper_nf"
THREADS=64

mkdir -p "$OUTDIR/bam" "$OUTDIR/logs"

module purge
module load bio/STAR/2.7.11b-GCC-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0
module load lang/Python/3.11.5-GCCcore-13.2.0

ulimit -n 65536

echo "════════════════════════════════════════════════════════"
echo "Mappability — v8_b_nf (MANE RNA reads, proper PE, 150 bp, dedup CDS (nf))"
echo "Reads from: $SIM_DIR"
echo "════════════════════════════════════════════════════════"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0d: Build Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna — ALL protein-coding CDS, genome-wide,
# deduplicated with seqkit (shared; skipped if already built)
#
# Split into two independently-resumable halves because seqkit has no
# confirmed module on this cluster: extraction (pure Python + pysam, no
# extra dependency) always runs here; dedup only runs here if seqkit is
# already on PATH, otherwise it stops with instructions to scp the raw
# FASTA down, run 'seqkit rmdup -s' locally, and scp Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna back
# to the same path here — the next run then just skips straight past this.
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0d: Checking Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna..."
if [ -f "$DEDUP_CDS_NF" ]; then
    echo "Step 0d: Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna already exists, skipping."
    SKIP_STEP_0D=1
else
    SKIP_STEP_0D=0
fi

if [ "$SKIP_STEP_0D" = "0" ] && [ -f "$DEDUP_CDS_NF_RAW" ]; then
    echo "Step 0d: Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_raw.fna already extracted, skipping extraction."
elif [ "$SKIP_STEP_0D" = "0" ]; then
    echo "Step 0d: Extracting genome-wide protein-coding CDS from nf GTF..."
python3 - << PYEOF
import pysam
from pathlib import Path
from collections import defaultdict

gtf_path  = Path("$GTF_NF")
genome_fa = Path("$GENOME_FA_NF")
raw_fa    = Path("$DEDUP_CDS_NF_RAW")

def parse_gtf_attrs(col9):
    attrs = defaultdict(list)
    for field in col9.strip().split(';'):
        field = field.strip()
        if not field:
            continue
        key, _, val = field.partition(' ')
        attrs[key].append(val.strip('"'))
    return attrs

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def revcomp(s): return s.translate(RC)[::-1]

print("  Scanning GTF for protein-coding CDS features...", flush=True)
tx_cds = defaultdict(list)
n_cds_lines = 0
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        a = parse_gtf_attrs(f[8])
        biotype = (a.get('transcript_biotype') or a.get('gene_biotype') or [''])[0]
        if biotype != 'protein_coding':
            continue
        tid = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None:
            continue
        versioned_tid = f"{tid}.{tver}" if tver else tid
        chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
        tx_cds[versioned_tid].append((chrom, start, end, strand))
        n_cds_lines += 1

for tid in tx_cds:
    tx_cds[tid].sort(key=lambda e: e[1])

print(f"  CDS feature lines: {n_cds_lines:,}")
print(f"  Protein-coding transcripts with CDS: {len(tx_cds):,}", flush=True)

print("  Fetching + splicing CDS sequences...", flush=True)
fa = pysam.FastaFile(str(genome_fa))
n_written = n_skip = 0
with open(raw_fa, 'w') as out:
    for tid, segs in tx_cds.items():
        strand = segs[0][3]
        try:
            seq = ''.join(fa.fetch(c, s, e) for (c, s, e, _st) in segs).upper()
        except (KeyError, ValueError):
            n_skip += 1
            continue
        if not seq:
            n_skip += 1
            continue
        if strand == '-':
            seq = revcomp(seq)
        out.write(f">{tid}\n{seq}\n")
        n_written += 1
fa.close()

print(f"  CDS sequences written (pre-dedup): {n_written:,}  (skipped: {n_skip})")
print(f"  Saved raw: {raw_fa}", flush=True)
PYEOF
fi

if [ "$SKIP_STEP_0D" = "0" ]; then
    if command -v seqkit &> /dev/null; then
        echo "  Deduplicating exact-sequence duplicates with seqkit rmdup -s..."
        seqkit rmdup -s -i "$DEDUP_CDS_NF_RAW" -o "$DEDUP_CDS_NF" --threads $THREADS
        echo "Step 0d complete."
    else
        echo ""
        echo "  seqkit not found on PATH on this cluster."
        echo "  $DEDUP_CDS_NF_RAW is ready (genome-wide protein-coding CDS, not yet deduplicated)."
        echo "  To finish Step 0d:"
        echo "    1. scp $DEDUP_CDS_NF_RAW to a machine with seqkit installed"
        echo "    2. run: seqkit rmdup -s -i Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_raw.fna -o Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna"
        echo "    3. scp Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna back to: $DEDUP_CDS_NF"
        echo "    4. re-run this script — it will detect the file and skip straight past Step 0d"
        exit 1
    fi
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0e: Build STAR index (dedup CDS, nf) — shared; skipped if already exists
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0e: Checking STAR dedup-CDS-nf index..."
mkdir -p "$CDS_IDX_NF"
if [ ! -f "$CDS_IDX_NF/SA" ]; then
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$CDS_IDX_NF" \
        --genomeFastaFiles    "$DEDUP_CDS_NF" \
        --genomeSAindexNbases 12 \
        --genomeChrBinNbits   11 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR/logs/star_index_cds_nf.log"
    echo "Step 0e complete."
else
    echo "Step 0e: index already exists, skipping."
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
# STEP 2: STAR alignment (PE, dedup CDS (nf))
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 2: STAR alignment (PE, dedup CDS (nf), $THREADS threads)..."

STAR \
    --runThreadN            $THREADS \
    --genomeDir             "$CDS_IDX_NF" \
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
out = outdir / "transcript_uniqueness_factors_dedup_cds_RNA_PE_L150bp_proper_nf.tsv"
fields = ['transcript_id','n_positions','n_unique_back','n_multi_back','n_unmapped','uniqueness_factor']
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(records)

fs = [r['uniqueness_factor'] for r in records]; n = len(fs)
print(f"  n={n:,}  mean={sum(fs)/n:.4f}  UF=1.0: {sum(1 for f in fs if f==1.0):,}  UF=0.0: {sum(1 for f in fs if f==0.0):,}")
print(f"  Saved: {out}")
PYEOF

echo "════════════════════════════════════════════════════════"
echo "Done: $OUTDIR/transcript_uniqueness_factors_dedup_cds_RNA_PE_L150bp_proper_nf.tsv"
echo "════════════════════════════════════════════════════════"
