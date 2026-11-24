#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — v2_b — NF REFERENCE (MANE-restricted dedup CDS)
#
# Reads     : exon-flank
# Reference : dedup CDS, MANE_Select-restricted (nf) — NOT the genome-wide
#             Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna built for v7_nf/v8_b_nf
# Mode      : SE
# Read len  : 100 bp
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

OUTDIR="$RESULTS/mappability_dedup_cds_SE_flank_nf"
THREADS=64
READ_LEN=100
FLANK=50

mkdir -p "$OUTDIR/simreads" "$OUTDIR/bam" "$OUTDIR/logs"

module purge
module load bio/STAR/2.7.11b-GCC-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0
module load lang/Python/3.11.5-GCCcore-13.2.0

ulimit -n 65536

echo "════════════════════════════════════════════════════════"
echo "Mappability — v2_b (nf reference, MANE-restricted dedup CDS, exon-flank, SE, 100 bp)"
echo "READ_LEN=$READ_LEN  THREADS=$THREADS"
echo "════════════════════════════════════════════════════════"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0: Expand chr_name_map.tsv -> chr_name_map_nf.tsv (does not touch original)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0: Expanding chr_name_map.tsv with Ensembl-style names..."
if [ -f "$RESULTS/chr_name_map_nf.tsv" ]; then
    echo "Step 0: chr_name_map_nf.tsv already exists, skipping."
else
if [ ! -f "$RESULTS/chr_name_map.tsv" ]; then
    echo "  $RESULTS/chr_name_map.tsv not found — rebuilding from the old NCBI genome FASTA"
    echo "  headers (same logic as v1_a.sh's own Step 0), instead of depending on a"
    echo "  pre-existing shared copy that no script actually guarantees exists."
    OLD_GENOME_FA_FOR_CHR_MAP="$REF/GCF_000001405.40_GRCh38.p14_genomic.fna"
python3 - << PYEOF
import re
from pathlib import Path

genome_fa = Path("$OLD_GENOME_FA_FOR_CHR_MAP")
out_map   = Path("$RESULTS/chr_name_map.tsv")

ncbi_map = {}
with open(genome_fa) as fh:
    for line in fh:
        if not line.startswith('>'):
            continue
        line = line.strip()
        name = line[1:].split()[0]
        rest = line[1:]
        m = re.search(r'chromosome (\w+)[,\s]', rest)
        if m and name.startswith('NC_'):
            ncbi_map[f"chr{m.group(1)}"] = name
        if 'mitochondrion' in rest.lower() or 'mitochondrial' in rest.lower():
            ncbi_map['chrM'] = name
            ncbi_map['chrMT'] = name

with open(out_map, 'w') as fh:
    for k, v in sorted(ncbi_map.items()):
        fh.write(f"{k}\t{v}\n")
print(f"  Chromosomes mapped: {len(ncbi_map)}")
print(f"  Wrote: {out_map}")
PYEOF
fi
python3 - << PYEOF
from pathlib import Path

old_map_path = Path("$RESULTS/chr_name_map.tsv")
genome_fa_nf = Path("$GENOME_FA_NF")
out_map      = Path("$RESULTS/chr_name_map_nf.tsv")

def ensembl_name_for(chr_name):
    if chr_name in ("chrM", "chrMT"):
        return "MT"
    if chr_name.startswith("chr"):
        return chr_name[3:]
    return chr_name

rows = []
with open(old_map_path) as fh:
    for line in fh:
        parts = line.strip().split('\t')
        if len(parts) != 2:
            continue
        chr_name, ncbi_name = parts
        rows.append([chr_name, ncbi_name, ensembl_name_for(chr_name), 'canonical'])

canonical_ensembl_names = {r[2] for r in rows}
seen = set()
with open(genome_fa_nf) as fh:
    for line in fh:
        if line.startswith('>'):
            seen.add(line[1:].split()[0])

for name in sorted(seen - canonical_ensembl_names):
    rows.append(['NA', 'NA', name, 'unplaced_scaffold_needs_crosswalk'])

rows.sort(key=lambda r: (r[3] != 'canonical', r[0]))
with open(out_map, 'w') as fh:
    fh.write("chr_name\tncbi_name\tensembl_name\tnote\n")
    for r in rows:
        fh.write('\t'.join(r) + '\n')
print(f"  Wrote: {out_map}")
PYEOF
fi
echo "Step 0 complete."

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
# STEP 1: Parse nf GTF (MANE_Select filter), fetch exon ± FLANK, simulate reads
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Simulating exon-flank reads (READ_LEN=$READ_LEN, FLANK=$FLANK, SE)..."

python3 - << PYEOF
import subprocess
from pathlib import Path
from collections import defaultdict

gtf_path  = Path("$GTF_NF")
genome_fa = Path("$GENOME_FA_NF")
outdir    = Path("$OUTDIR")
read_len  = int("$READ_LEN")
flank     = int("$FLANK")


def parse_gtf_attrs(col9):
    attrs = defaultdict(list)
    for field in col9.strip().split(';'):
        field = field.strip()
        if not field: continue
        key, _, val = field.partition(' ')
        attrs[key].append(val.strip('"'))
    return attrs

mane_select_tx = set()
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'transcript': continue
        a = parse_gtf_attrs(f[8])
        if 'MANE_Select' not in a.get('tag', []): continue
        tid = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None: continue
        mane_select_tx.add(f"{tid}.{tver}" if tver else tid)
print(f"  MANE_Select transcripts: {len(mane_select_tx):,}", flush=True)

tx_exons = defaultdict(list)
n_exons = 0
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        a = parse_gtf_attrs(f[8])
        tid = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None: continue
        versioned_tid = f"{tid}.{tver}" if tver else tid
        if versioned_tid not in mane_select_tx: continue
        chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
        tx_exons[versioned_tid].append((chrom, start, end, strand))
        n_exons += 1
print(f"  Exon features parsed: {n_exons:,}", flush=True)

fai = Path(str(genome_fa) + '.fai')
if not fai.exists():
    subprocess.run(['samtools', 'faidx', str(genome_fa)], check=True)

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
# STEP 2: STAR alignment — SE, MANE-restricted dedup CDS (nf)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 2: STAR alignment (SE, dedup-CDS-MANE-nf, $THREADS threads)..."

STAR \
    --runThreadN            $THREADS \
    --genomeDir             "$CDS_MANE_IDX_NF" \
    --readFilesIn           "$OUTDIR/simreads/sim_R1.fastq" \
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
    if read.is_supplementary or read.is_secondary or read.is_unmapped: continue
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
out = outdir / "transcript_uniqueness_factors_dedup_cds_SE_flank_L100bp_nf.tsv"
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
echo "Done: $OUTDIR/transcript_uniqueness_factors_dedup_cds_SE_flank_L100bp_nf.tsv"
echo "Compare against: results/transcript_uniqueness_factors_dedup_cds_SE_flank_L100bp.tsv"
echo "════════════════════════════════════════════════════════"
