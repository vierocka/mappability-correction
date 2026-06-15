#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=kallisto_crossval_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL

# ─────────────────────────────────────────────────────────────────────────────
# Kallisto cross-validation against MANE RNA and deduplicated CDS.
#
# Three read sources:
#   MANE spliced RNA reads  — reused from STAR v5 scripts (polyA proxy)
#   Exon-flank reads        — reused from STAR v5 scripts; intronic content;
#                             intron retention proxy; expected to fail pseudoalignment
#   Dedup CDS reads         — simulated here from the dedup CDS FASTA (93,088 seqs);
#                             quantifies CDS-level sequence uniqueness and coding
#                             redundancy across the annotated genome independently
#                             of UTR content
#
# Two Kallisto indices:
#   MANE RNA   — 19,437 full spliced transcripts (UTRs included)
#   Dedup CDS  — 93,088 CDS sequences (no UTRs; sequence-unique pseudogenes retained)
#
# Full genome Kallisto index excluded — not designed for genome-scale references.
# No bootstraps (-b 0) — not needed for mappability assessment.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$SCRIPT_DIR/../../Ref"
RESULTS="$SCRIPT_DIR/results"
THREADS=24

MANE_FA="$REF/MANE.GRCh38.v1.5.ensembl_rna.fna"
CDS_FA="$REF/GCF_000001405.40_GRCh38.p14_cds_from_genomic.dedup.fna"
MANE_IDX="$REF/kallisto_index_MANE.idx"
CDS_IDX="$REF/kallisto_index_dedup_cds.idx"

echo "════════════════════════════════════════════════════════"
echo "Kallisto cross-validation — mappability"
echo "Started : $(date)"
echo "Host    : $(hostname)"
echo "════════════════════════════════════════════════════════"

module purge
module load bio/kallisto/0.50.1-gompi-2022b

# ── Step 0: Build indices ─────────────────────────────────────────────────────

echo ""
echo "Step 0: Building Kallisto indices..."

if [ ! -f "$MANE_IDX" ]; then
    echo "  Building MANE RNA index..."
    kallisto index -i "$MANE_IDX" "$MANE_FA"
    echo "  Done."
else
    echo "  MANE RNA index exists, skipping."
fi

if [ ! -f "$CDS_IDX" ]; then
    echo "  Building deduplicated CDS index..."
    kallisto index -i "$CDS_IDX" "$CDS_FA"
    echo "  Done."
else
    echo "  Deduplicated CDS index exists, skipping."
fi

# ── Step 1: Simulate reads from deduplicated CDS ──────────────────────────────
#
# Reads from MANE and exon-flank sources already exist (produced by STAR v5
# scripts). Dedup CDS reads do not exist elsewhere and are produced here.
# SE and PE at four read lengths: 75, 100, 150, 200 bp.

echo ""
echo "Step 1: Simulating reads from deduplicated CDS FASTA..."

simulate_cds() {
    local read_len="$1"
    local mode="$2"     # SE or PE
    local outdir="$3"

    mkdir -p "$outdir/simreads"

    if [ -f "$outdir/simreads/sim_R1.fastq" ]; then
        echo "  Reads exist, skipping ($outdir)."
        return
    fi

    echo "  Simulating $mode ${read_len} bp → $(basename "$outdir")..."
    python3 - "$CDS_FA" "$outdir/simreads" "$read_len" "$mode" << 'PYEOF'
import sys
from pathlib import Path

cds_fa = Path(sys.argv[1])
simdir = Path(sys.argv[2])
rl     = int(sys.argv[3])
mode   = sys.argv[4]

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def revcomp(s): return s.translate(RC)[::-1]

def parse_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name: yield name, ''.join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.upper())
    if name: yield name, ''.join(seq)

qual    = 'I' * rl
n_reads = 0
n_skip  = 0
n_seq   = 0

r1_path = simdir / 'sim_R1.fastq'
r2_path = simdir / 'sim_R2.fastq'

r2 = open(r2_path, 'w') if mode == 'PE' else None

with open(r1_path, 'w') as r1:
    for seq_id, seq in parse_fasta(cds_fa):
        L = len(seq)
        if L < rl:
            n_skip += 1
            continue
        n_seq += 1
        for pos in range(L - rl + 1):
            fwd = seq[pos : pos + rl]
            if fwd.count('N') > rl // 5:
                continue
            name = f"{seq_id}|{pos}"
            r1.write(f"@{name}\n{fwd}\n+\n{qual}\n")
            if r2:
                r2.write(f"@{name}\n{revcomp(fwd)}\n+\n{qual}\n")
            n_reads += 1

if r2:
    r2.close()

print(f"  Sequences: {n_seq:,}  (too short: {n_skip})")
print(f"  Reads simulated: {n_reads:,}", flush=True)
PYEOF
}

for LEN in 75 100 150 200; do
    SUFFIX=""
    [ "$LEN" != "100" ] && SUFFIX="_L${LEN}bp"

    simulate_cds "$LEN" SE "$RESULTS/mappability_dedup_cds_reads_SE${SUFFIX}"
    simulate_cds "$LEN" PE "$RESULTS/mappability_dedup_cds_reads_PE${SUFFIX}"
done

echo "Step 1 complete."

# ── Helper functions for kallisto quant ───────────────────────────────────────

run_se() {
    local outdir="$1"
    local read_len="$2"
    local r1="$outdir/simreads/sim_R1.fastq"

    if [ ! -f "$r1" ]; then
        echo "  WARNING: reads not found — $r1 — skipping."
        return
    fi

    echo "  → MANE RNA index"
    mkdir -p "$outdir/kallisto_vs_MANE"
    kallisto quant \
        --single -l "$read_len" -s 1 \
        -b 0 -t "$THREADS" \
        -i "$MANE_IDX" \
        -o "$outdir/kallisto_vs_MANE" \
        "$r1"

    echo "  → Dedup CDS index"
    mkdir -p "$outdir/kallisto_vs_dedup_cds"
    kallisto quant \
        --single -l "$read_len" -s 1 \
        -b 0 -t "$THREADS" \
        -i "$CDS_IDX" \
        -o "$outdir/kallisto_vs_dedup_cds" \
        "$r1"
}

run_pe() {
    local outdir="$1"
    local r1="$outdir/simreads/sim_R1.fastq"
    local r2="$outdir/simreads/sim_R2.fastq"

    if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
        echo "  WARNING: reads not found — $r1 / $r2 — skipping."
        return
    fi

    echo "  → MANE RNA index"
    mkdir -p "$outdir/kallisto_vs_MANE"
    kallisto quant \
        -b 0 -t "$THREADS" \
        -i "$MANE_IDX" \
        -o "$outdir/kallisto_vs_MANE" \
        "$r1" "$r2"

    echo "  → Dedup CDS index"
    mkdir -p "$outdir/kallisto_vs_dedup_cds"
    kallisto quant \
        -b 0 -t "$THREADS" \
        -i "$CDS_IDX" \
        -o "$outdir/kallisto_vs_dedup_cds" \
        "$r1" "$r2"
}

# ── MANE spliced RNA reads — SE ───────────────────────────────────────────────

echo ""
echo "━━━━ MANE RNA reads — SE ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "v5_c :  75 bp"
run_se "$RESULTS/mappability_MANE_RNA_SE_L75bp"  75

echo "v5_a : 100 bp"
run_se "$RESULTS/mappability_MANE_RNA_SE"        100

echo "v5_e : 150 bp"
run_se "$RESULTS/mappability_MANE_RNA_SE_L150bp" 150

echo "v5_d : 200 bp"
run_se "$RESULTS/mappability_MANE_RNA_SE_L200bp" 200

# ── MANE spliced RNA reads — PE ───────────────────────────────────────────────

echo ""
echo "━━━━ MANE RNA reads — PE ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "v5_f :  75 bp"
run_pe "$RESULTS/mappability_MANE_RNA_PE_L75bp"

echo "v5_g : 100 bp"
run_pe "$RESULTS/mappability_MANE_RNA_PE"

echo "v5_h : 150 bp"
run_pe "$RESULTS/mappability_MANE_RNA_PE_L150bp"

echo "v5_i : 200 bp"
run_pe "$RESULTS/mappability_MANE_RNA_PE_L200bp"

# ── Exon-flank reads — SE (intron retention proxy) ───────────────────────────

echo ""
echo "━━━━ Exon-flank reads — SE (intron retention proxy) ━━━━━━━━━━━━━━━━━━━━━"

echo "v5_j :  75 bp"
run_se "$RESULTS/mappability_MANE_exonflank_SE_L75bp"  75

echo "v5_k : 100 bp"
run_se "$RESULTS/mappability_MANE_exonflank_SE"        100

echo "v5_l : 150 bp"
run_se "$RESULTS/mappability_MANE_exonflank_SE_L150bp" 150

# ── Exon-flank reads — PE (intron retention proxy) ───────────────────────────

echo ""
echo "━━━━ Exon-flank reads — PE (intron retention proxy) ━━━━━━━━━━━━━━━━━━━━━"

echo "v5_m :  75 bp"
run_pe "$RESULTS/mappability_MANE_exonflank_PE_L75bp"

echo "v5_b : 100 bp"
run_pe "$RESULTS/mappability_MANE_exonflank_PE"

echo "v5_n : 150 bp"
run_pe "$RESULTS/mappability_MANE_exonflank_PE_L150bp"

# ── Dedup CDS reads — SE ─────────────────────────────────────────────────────

echo ""
echo "━━━━ Dedup CDS reads — SE (CDS uniqueness / coding redundancy) ━━━━━━━━━━"

echo "75 bp"
run_se "$RESULTS/mappability_dedup_cds_reads_SE_L75bp"  75

echo "100 bp"
run_se "$RESULTS/mappability_dedup_cds_reads_SE"        100

echo "150 bp"
run_se "$RESULTS/mappability_dedup_cds_reads_SE_L150bp" 150

echo "200 bp"
run_se "$RESULTS/mappability_dedup_cds_reads_SE_L200bp" 200

# ── Dedup CDS reads — PE ─────────────────────────────────────────────────────

echo ""
echo "━━━━ Dedup CDS reads — PE (CDS uniqueness / coding redundancy) ━━━━━━━━━━"

echo "75 bp"
run_pe "$RESULTS/mappability_dedup_cds_reads_PE_L75bp"

echo "100 bp"
run_pe "$RESULTS/mappability_dedup_cds_reads_PE"

echo "150 bp"
run_pe "$RESULTS/mappability_dedup_cds_reads_PE_L150bp"

echo "200 bp"
run_pe "$RESULTS/mappability_dedup_cds_reads_PE_L200bp"

# ── Summary ───────────────────────────────────────────────────────────────────

echo ""
echo "════════════════════════════════════════════════════════"
echo "All Kallisto runs complete: $(date)"
echo ""
echo "Result files (abundance.tsv):"
find "$RESULTS" -name "abundance.tsv" -path "*/kallisto_vs_*" | sort
echo "════════════════════════════════════════════════════════"
