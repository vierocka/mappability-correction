#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — genomic exon-flank approach — NF REFERENCE
#
# This is the "_nf" counterpart of mappability_correction_v1_a.sh: same method
# (exon ± flank windows, PE, 100 bp), but aligned against the nf-core-style
# Ensembl reference instead of the NCBI RefSeq GRCh38.p14 reference, so the two
# can be compared transcript-by-transcript.
#
# Reference source (Ensembl, GRCh38.116):
#   Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa — numbered chromosomes
#   (1,2,...,22,X,Y,MT) + unplaced/unlocalized scaffolds, no alt-loci haplotypes.
#   Homo_sapiens.GRCh38.116.gtf — exon coordinates for the SAME assembly,
#   already using the identical chromosome names as the genome FASTA.
#
# Transcript scope: kept identical to the rest of the pipeline by filtering the
# GTF to transcripts carrying `tag "MANE_Select"` (~19.4k transcripts, same set
# used by MANE.GRCh38.v1.5 elsewhere). Ensembl GTF repeats `tag` as multiple
# attributes on one line (e.g. tag "basic"; tag "MANE_Select";), so parsing
# collects ALL tag values per line rather than taking the first one.
#
# Because the GTF and genome FASTA share chromosome names natively, no chr→NC
# translation is needed here (unlike v1_a.sh, which needs chr_name_map.tsv to
# go from MANE GFF's "chr1" to NCBI's "NC_000001.11"). Step 0 below still
# EXPANDS results/chr_name_map.tsv into results/chr_name_map_nf.tsv (adding an
# ensembl_name column) — that expanded map is not used internally by this
# script; it exists purely so results from the old and new reference can be
# cross-referenced by chromosome later.
#
# Nothing from the existing pipeline is modified or reused destructively:
#   - new STAR index   : $REF/STAR_index_genome_nf       (old STAR_index_genome untouched)
#   - new output dir   : results/mappability_genomic_nf  (old mappability_genomic untouched)
#   - new chr map      : results/chr_name_map_nf.tsv      (old chr_name_map.tsv untouched)
#   - new result TSV   : transcript_uniqueness_factors_genomic_L100bp_nf.tsv
#
# Before running this at scale, run (from the mappability/ dir):
#   python3 check_id_conversion_nf.py   --gtf <path-to-Ref-GTF>
#   python3 check_exon_positions_nf.py  --new-gtf <path-to-Ref-GTF>
# to confirm the MANE_Select transcript set and exon coordinates agree with
# the existing MANE-based pipeline before trusting a mapping-difference result
# as being due to the reference and not an annotation mismatch.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS="$SCRIPT_DIR/../results"

# ── nf reference (Ensembl GRCh38.116) — lives in Ref/ alongside the original
#    NCBI RefSeq files, distinguished by filename, not a separate directory ──
REF="$SCRIPT_DIR/../../Ref"
GENOME_FA_NF="$REF/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
GTF_NF="$REF/Homo_sapiens.GRCh38.116.gtf"
GENOME_IDX_NF="$REF/STAR_index_genome_nf"

OUTDIR="$RESULTS/mappability_genomic_nf"
THREADS=64
READ_LEN=100
FLANK=50

mkdir -p "$OUTDIR/simreads" "$OUTDIR/bam" "$OUTDIR/logs"

# ── Load modules (standalone script — same set as run_all_mappability.sh /
#    mappability_v8_HPC.sh / run_3prime_HPC.sh) ───────────────────────────────
module purge
module load bio/STAR/2.7.11b-GCC-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0
module load lang/Python/3.11.5-GCCcore-13.2.0

# STAR BAM sorting needs many open file descriptors
ulimit -n 65536

echo "════════════════════════════════════════════════════════"
echo "MANE mappability — genomic exon-flank approach — nf reference"
echo "READ_LEN=$READ_LEN  FLANK=$FLANK  THREADS=$THREADS"
echo "════════════════════════════════════════════════════════"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0: Expand chr_name_map.tsv → chr_name_map_nf.tsv (does not touch original)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0: Expanding chr_name_map.tsv with Ensembl-style names..."

# results/chr_name_map.tsv is NOT produced by any script — v1_a.sh writes its
# chr->NC map to its own per-run OUTDIR, never to a shared results/ location.
# Rather than depend on someone having copied a shared one here by hand,
# rebuild it from the old NCBI genome FASTA headers if it's missing, using
# the same header-parsing logic as v1_a.sh's own Step 0.
if [ ! -f "$RESULTS/chr_name_map.tsv" ]; then
    echo "  $RESULTS/chr_name_map.tsv not found — rebuilding from the old NCBI genome FASTA"
    echo "  headers (same logic as v1_a.sh's own Step 0)."
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
out_map_run  = Path("$OUTDIR/chr_name_map_nf.tsv")

def ensembl_name_for(chr_name):
    if chr_name in ("chrM", "chrMT"):
        return "MT"
    if chr_name.startswith("chr"):
        return chr_name[3:]
    return chr_name

rows = []  # (chr_name, ncbi_name, ensembl_name, note)
with open(old_map_path) as fh:
    for line in fh:
        parts = line.strip().split('\t')
        if len(parts) != 2:
            continue
        chr_name, ncbi_name = parts
        rows.append([chr_name, ncbi_name, ensembl_name_for(chr_name), 'canonical'])

canonical_ensembl_names = {r[2] for r in rows}

print("  Scanning nf genome FASTA headers for additional contigs...", flush=True)
seen = set()
with open(genome_fa_nf) as fh:
    for line in fh:
        if not line.startswith('>'):
            continue
        name = line[1:].split()[0]
        seen.add(name)

extra = sorted(seen - canonical_ensembl_names)
for name in extra:
    rows.append(['NA', 'NA', name, 'unplaced_scaffold_needs_crosswalk'])

rows.sort(key=lambda r: (r[3] != 'canonical', r[0]))

for path in (out_map, out_map_run):
    with open(path, 'w') as fh:
        fh.write("chr_name\tncbi_name\tensembl_name\tnote\n")
        for r in rows:
            fh.write('\t'.join(r) + '\n')

n_canon = sum(1 for r in rows if r[3] == 'canonical')
n_extra = len(rows) - n_canon
print(f"  Canonical chromosomes mapped : {n_canon}")
print(f"  Unmatched scaffolds found    : {n_extra}"
      f"{' (see note column — no RefSeq accession crosswalk attempted)' if n_extra else ''}")
print(f"  Wrote: {out_map}")
print(f"  Wrote: {out_map_run}")
PYEOF
echo "Step 0 complete."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0b: Build a NEW STAR index from the nf genome (never touches the old one)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0b: Building STAR index from nf genome (new index, old index untouched)..."
mkdir -p "$GENOME_IDX_NF"
if [ ! -f "$GENOME_IDX_NF/SA" ]; then
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$GENOME_IDX_NF" \
        --genomeFastaFiles    "$GENOME_FA_NF" \
        --genomeSAindexNbases 14 \
        --genomeChrBinNbits   14 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR/logs/star_index_nf.log"
    echo "Step 0b complete."
else
    echo "Step 0b: nf index already exists, skipping."
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Parse GTF (MANE_Select-tagged transcripts only), fetch exon ± FLANK,
#         simulate reads
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Parsing nf GTF (MANE_Select filter) and simulating reads..."

python3 - << PYEOF
import subprocess
from pathlib import Path
from collections import defaultdict

gtf_path  = Path("$GTF_NF")
genome_fa = Path("$GENOME_FA_NF")
outdir    = Path("$OUTDIR")
read_len  = $READ_LEN
flank     = $FLANK

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def revcomp(s): return s.translate(RC)[::-1]

def parse_gtf_attrs(col9):
    """Ensembl GTF attributes: key "value"; repeated keys (esp. 'tag') accumulate."""
    attrs = defaultdict(list)
    for field in col9.strip().split(';'):
        field = field.strip()
        if not field:
            continue
        key, _, val = field.partition(' ')
        attrs[key].append(val.strip('"'))
    return attrs

# ── Pass 1: find MANE_Select-tagged transcripts ───────────────────────────────
print("  Pass 1: scanning for MANE_Select-tagged transcripts...", flush=True)
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

print(f"  MANE_Select transcripts found: {len(mane_select_tx):,}")

# ── Pass 2: collect exon coordinates for those transcripts ────────────────────
print("  Pass 2: collecting exon coordinates...", flush=True)
tx_exons = defaultdict(list)
n_exons = 0
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        a = parse_gtf_attrs(f[8])
        tid = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None:
            continue
        versioned_tid = f"{tid}.{tver}" if tver else tid
        if versioned_tid not in mane_select_tx:
            continue
        chrom  = f[0]                 # already matches genome FASTA names natively
        start  = int(f[3]) - 1        # GTF 1-based -> 0-based
        end    = int(f[4])
        strand = f[6]
        tx_exons[versioned_tid].append((chrom, start, end, strand))
        n_exons += 1

print(f"  Exon features parsed: {n_exons:,}")
print(f"  Transcripts with exons: {len(tx_exons):,}", flush=True)

# ── Index genome FASTA ────────────────────────────────────────────────────────
fai = Path(str(genome_fa) + '.fai')
if not fai.exists():
    print("  Indexing genome with samtools faidx...", flush=True)
    subprocess.run(['samtools', 'faidx', str(genome_fa)], check=True)

# ── Simulate reads from exon ± flank windows ─────────────────────────────────
print(f"  Simulating reads (step=1, READ_LEN={read_len}, FLANK={flank})...", flush=True)

r1_path = outdir / "simreads/sim_R1.fastq"
r2_path = outdir / "simreads/sim_R2.fastq"
qual    = 'I' * read_len

n_reads = 0
n_skip  = 0

with open(r1_path, 'w') as r1, open(r2_path, 'w') as r2:
    for tid, exons in tx_exons.items():
        for ex_idx, (chrom, ex_start, ex_end, strand) in enumerate(exons):
            win_start = max(0, ex_start - flank)
            win_end   = ex_end + flank
            region = f"{chrom}:{win_start + 1}-{win_end}"
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
# STEP 2: STAR alignment — splice-aware (introns enabled) — against nf index
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 2: STAR alignment against nf genome index ($THREADS threads)..."

STAR \
    --runThreadN      $THREADS \
    --genomeDir       "$GENOME_IDX_NF" \
    --readFilesIn     "$OUTDIR/simreads/sim_R1.fastq" \
                      "$OUTDIR/simreads/sim_R2.fastq" \
    --outSAMtype      BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM \
    --outSAMmultNmax  1 \
    --outFilterMultimapNmax 20 \
    --outBAMsortingThreadN $THREADS \
    --outBAMsortingBinsN   20 \
    --limitBAMsortRAM      160000000000 \
    --outFileNamePrefix "$OUTDIR/bam/sim_" \
    2>&1 | tee "$OUTDIR/logs/star_sim_nf.log"

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

print("  Counting simulated reads per transcript...", flush=True)
sim_counts = defaultdict(int)
with open(r1_path) as fh:
    for line in fh:
        if line.startswith('@'):
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
    if read.is_unmapped:
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

records.sort(key=lambda r: r['uniqueness_factor'])
out = outdir / "transcript_uniqueness_factors_genomic_L100bp_nf.tsv"
fields = ['transcript_id','n_positions','n_unique_back','n_multi_back','n_unmapped','uniqueness_factor']
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(records)

fs = [r['uniqueness_factor'] for r in records]
n  = len(fs)
print(f"\n  ── Distribution ─────────────────────────────────────")
print(f"  n         : {n:,}")
print(f"  mean      : {sum(fs)/n:.4f}" if n else "  mean      : n/a")
print(f"  f=1.00    : {sum(1 for x in fs if x==1.0):,}")
print(f"  f>=0.90   : {sum(1 for x in fs if x>=0.90):,}")
print(f"  f 0.10-0.90: {sum(1 for x in fs if 0.10<=x<1.0):,}")
print(f"  f<0.10    : {sum(1 for x in fs if x<0.10):,}")
print(f"  f=0.00    : {sum(1 for x in fs if x==0.0):,}")
print(f"\n  Saved: {out}")
PYEOF

echo ""
echo "════════════════════════════════════════════════════════"
echo "Done: $OUTDIR/transcript_uniqueness_factors_genomic_L100bp_nf.tsv"
echo "Compare against: results/mappability_genomic/transcript_uniqueness_factors_genomic_L100bp.tsv"
echo "════════════════════════════════════════════════════════"
