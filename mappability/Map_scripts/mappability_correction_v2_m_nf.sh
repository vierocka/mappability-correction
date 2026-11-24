#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — v2_m — NF REFERENCE (MANE-restricted dedup CDS)
#
# Reads     : MANE RNA
# Reference : dedup CDS, MANE_Select-restricted (nf) — NOT the genome-wide
#             Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna built for v7_nf/v8_b_nf
# Mode      : PE
# Read len  : 150 bp
# alignIntronMax=1 : False
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS="$SCRIPT_DIR/../results"

REF="$SCRIPT_DIR/../../Ref"
GENOME_FA_NF="$REF/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
GTF_NF="$REF/Homo_sapiens.GRCh38.116.gtf"
MANE_RNA_NF="$REF/MANE_RNA_nf.fna"
DEDUP_CDS_MANE_NF_RAW="$REF/Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane_raw.fna"
DEDUP_CDS_MANE_NF="$REF/Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna"
CDS_MANE_IDX_NF="$REF/STAR_index_dedup_cds_mane_nf"

OUTDIR="$RESULTS/mappability_dedup_cds_RNA_PE_L150bp_nf"
THREADS=64
READ_LEN=150
FLANK=50

mkdir -p "$OUTDIR/simreads" "$OUTDIR/bam" "$OUTDIR/logs"

module purge
module load bio/STAR/2.7.11b-GCC-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0
module load lang/Python/3.11.5-GCCcore-13.2.0

ulimit -n 65536

echo "════════════════════════════════════════════════════════"
echo "Mappability — v2_m (nf reference, MANE-restricted dedup CDS, MANE RNA, PE, 150 bp)"
echo "READ_LEN=$READ_LEN  THREADS=$THREADS"
echo "════════════════════════════════════════════════════════"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0d: Build Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna — MANE_Select-restricted CDS,
# deduplicated with seqkit (shared; skipped if already built)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0d: Checking Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna..."
if [ -f "$DEDUP_CDS_MANE_NF" ]; then
    echo "Step 0d: Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna already exists, skipping."
    SKIP_STEP_0D=1
else
    SKIP_STEP_0D=0
fi

if [ "$SKIP_STEP_0D" = "0" ] && [ -f "$DEDUP_CDS_MANE_NF_RAW" ]; then
    echo "Step 0d: Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane_raw.fna already extracted, skipping extraction."
elif [ "$SKIP_STEP_0D" = "0" ]; then
    echo "Step 0d: Extracting MANE_Select CDS from nf GTF..."
python3 - << PYEOF
import pysam
from pathlib import Path
from collections import defaultdict

gtf_path  = Path("$GTF_NF")
genome_fa = Path("$GENOME_FA_NF")
raw_fa    = Path("$DEDUP_CDS_MANE_NF_RAW")

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

mane_select_tx = set()
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'transcript':
            continue
        a = parse_gtf_attrs(f[8])
        if 'MANE_Select' not in a.get('tag', []):
            continue
        tid = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None:
            continue
        mane_select_tx.add(f"{tid}.{tver}" if tver else tid)
print(f"  MANE_Select transcripts: {len(mane_select_tx):,}", flush=True)

tx_cds = defaultdict(list)
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        a = parse_gtf_attrs(f[8])
        tid = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None:
            continue
        versioned_tid = f"{tid}.{tver}" if tver else tid
        if versioned_tid not in mane_select_tx:
            continue
        chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
        tx_cds[versioned_tid].append((chrom, start, end, strand))

for tid in tx_cds:
    tx_cds[tid].sort(key=lambda e: e[1])
print(f"  MANE_Select transcripts with CDS: {len(tx_cds):,}", flush=True)

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
        seqkit rmdup -s -i "$DEDUP_CDS_MANE_NF_RAW" -o "$DEDUP_CDS_MANE_NF" --threads $THREADS
        echo "Step 0d complete."
    else
        echo ""
        echo "  seqkit not found on PATH on this cluster."
        echo "  $DEDUP_CDS_MANE_NF_RAW is ready (MANE_Select CDS, not yet deduplicated)."
        echo "  To finish Step 0d:"
        echo "    1. scp $DEDUP_CDS_MANE_NF_RAW to a machine with seqkit installed"
        echo "    2. run: seqkit rmdup -s -i Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane_raw.fna -o Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna"
        echo "    3. scp Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna back to: $DEDUP_CDS_MANE_NF"
        echo "    4. re-run this script — it will detect the file and skip straight past Step 0d"
        exit 1
    fi
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0e: Build STAR index (dedup CDS, MANE-restricted, nf) — shared; skipped if exists
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0e: Checking STAR dedup-CDS-MANE-nf index..."
mkdir -p "$CDS_MANE_IDX_NF"
if [ ! -f "$CDS_MANE_IDX_NF/SA" ]; then
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$CDS_MANE_IDX_NF" \
        --genomeFastaFiles    "$DEDUP_CDS_MANE_NF" \
        --genomeSAindexNbases 12 \
        --genomeChrBinNbits   11 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR/logs/star_index_cds_mane_nf.log"
    echo "Step 0e complete."
else
    echo "Step 0e: index already exists, skipping."
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0c: Build MANE_RNA_nf.fna — spliced MANE_Select transcript sequences,
# derived from the nf genome + GTF (nf-core recommends deriving this rather
# than downloading a separate transcriptome FASTA). Shared; skipped if exists.
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0c: Checking MANE_RNA_nf.fna..."
if [ -f "$MANE_RNA_NF" ]; then
    echo "Step 0c: MANE_RNA_nf.fna already exists, skipping."
else
    echo "Step 0c: Splicing MANE_Select transcripts from nf GTF..."
    if [ ! -f "$GENOME_FA_NF.fai" ]; then
        samtools faidx "$GENOME_FA_NF"
    fi
python3 - << PYEOF
import subprocess
from pathlib import Path
from collections import defaultdict

gtf_path  = Path("$GTF_NF")
genome_fa = Path("$GENOME_FA_NF")
out_fa    = Path("$MANE_RNA_NF")

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def revcomp(s): return s.translate(RC)[::-1]

def parse_gtf_attrs(col9):
    attrs = defaultdict(list)
    for field in col9.strip().split(';'):
        field = field.strip()
        if not field:
            continue
        key, _, val = field.partition(' ')
        attrs[key].append(val.strip('"'))
    return attrs

mane_select_tx = set()
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'transcript':
            continue
        a = parse_gtf_attrs(f[8])
        if 'MANE_Select' not in a.get('tag', []):
            continue
        tid  = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None:
            continue
        mane_select_tx.add(f"{tid}.{tver}" if tver else tid)
print(f"  MANE_Select transcripts: {len(mane_select_tx):,}", flush=True)

tx_exons = defaultdict(list)
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        a = parse_gtf_attrs(f[8])
        tid  = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None:
            continue
        versioned_tid = f"{tid}.{tver}" if tver else tid
        if versioned_tid not in mane_select_tx:
            continue
        tx_exons[versioned_tid].append((f[0], int(f[3]) - 1, int(f[4]), f[6]))

for tid in tx_exons:
    tx_exons[tid].sort(key=lambda e: e[1])

n_written = n_skip = 0
with open(out_fa, 'w') as out:
    for tid, exons in tx_exons.items():
        chrom  = exons[0][0]
        strand = exons[0][3]
        seq_parts = []
        ok = True
        for (c, start, end, _s) in exons:
            region = f"{c}:{start + 1}-{end}"
            try:
                res = subprocess.run(['samtools', 'faidx', str(genome_fa), region],
                                     capture_output=True, text=True, check=True)
                seq_parts.append(''.join(res.stdout.split('\n')[1:]).upper())
            except subprocess.CalledProcessError:
                ok = False
                break
        if not ok:
            n_skip += 1
            continue
        seq = ''.join(seq_parts)
        if strand == '-':
            seq = revcomp(seq)
        out.write(f">{tid}\n{seq}\n")
        n_written += 1

print(f"  Transcripts written: {n_written:,}  (skipped: {n_skip})")
print(f"  Saved: {out_fa}", flush=True)
PYEOF
fi
echo "Step 0c complete."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Simulate reads from spliced MANE_RNA_nf.fna
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Simulating MANE RNA reads (READ_LEN=$READ_LEN, PE)..."

python3 - << PYEOF
from pathlib import Path

mane_rna = Path("$MANE_RNA_NF")
outdir   = Path("$OUTDIR")
read_len = int("$READ_LEN")

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

r1_path = outdir / 'simreads/sim_R1.fastq'
r2_path = outdir / 'simreads/sim_R2.fastq'
qual = 'I' * read_len
n_reads = n_tx = n_skip = 0

with open(r1_path, 'w') as r1, open(r2_path, 'w') as r2:
    for tid, seq in parse_fasta(mane_rna):
        if len(seq) < read_len: n_skip += 1; continue
        n_tx += 1
        for pos in range(len(seq) - read_len + 1):
            fwd = seq[pos : pos + read_len]
            if fwd.count('N') > read_len // 5: continue
            r1.write(f'@{tid}|{pos}\n{fwd}\n+\n{qual}\n')
            r2.write(f'@{tid}|{pos}\n{revcomp(fwd)}\n+\n{qual}\n')
            n_reads += 1

print(f"  Transcripts: {n_tx:,}  (skipped too short: {n_skip})")
print(f"  Total reads simulated: {n_reads:,}", flush=True)
PYEOF
echo "Step 1 complete."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: STAR alignment — PE, MANE-restricted dedup CDS (nf)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 2: STAR alignment (PE, dedup-CDS-MANE-nf, $THREADS threads)..."

STAR \
    --runThreadN            $THREADS \
    --genomeDir             "$CDS_MANE_IDX_NF" \
    --readFilesIn           "$OUTDIR/simreads/sim_R1.fastq" \
                            "$OUTDIR/simreads/sim_R2.fastq" \
    --outSAMtype            BAM SortedByCoordinate \
    --outSAMattributes      NH HI AS NM \
    --outSAMmultNmax        1 \
    --outFilterMultimapNmax 20 \
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
print(f"  Transcripts: {len(sim_counts):,}  Reads: {sum(sim_counts.values()):,}", flush=True)

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
out = outdir / "transcript_uniqueness_factors_dedup_cds_RNA_PE_L150bp_nf.tsv"
fields = ['transcript_id','n_positions','n_unique_back','n_multi_back','n_unmapped','uniqueness_factor']
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(records)

fs = [r['uniqueness_factor'] for r in records]; n = len(fs)
print(f"\n  ── Distribution ─────────────────────────────────────")
print(f"  n         : {n:,}")
if n:
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
echo "Done: $OUTDIR/transcript_uniqueness_factors_dedup_cds_RNA_PE_L150bp_nf.tsv"
echo "Compare against: results/transcript_uniqueness_factors_dedup_cds_RNA_PE_L150bp.tsv"
echo "════════════════════════════════════════════════════════"
