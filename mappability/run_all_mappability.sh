#!/bin/bash -l
#SBATCH --cpus-per-task=64
#SBATCH --mem=192gb
#SBATCH --time=72:00:00
#SBATCH --output=mappability_all_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=vkovacov@uni-koeln.de
#SBATCH --account=ag-laessig

# ─────────────────────────────────────────────────────────────────────────────
# Run all 9 mappability correction scripts sequentially.
#
# Execution order matters:
#   v1_a runs first  → builds the shared genome STAR index (STAR_index_genome)
#   v2_a runs first of the CDS group → builds STAR_index_dedup_cds
#   All subsequent scripts reuse whichever index already exists.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

echo "════════════════════════════════════════════════════════"
echo "Mappability correction — full run (13 scripts)"
echo "Started : $(date)"
echo "Host    : $(hostname)"
echo "Dir     : $SCRIPT_DIR"
echo "════════════════════════════════════════════════════════"

# ── Load modules ──────────────────────────────────────────────────────────────
module purge
module load bio/STAR/2.7.11b-GCC-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0
module load lang/Python/3.11.5-GCCcore-13.2.0

# ── Raise open-file limit (STAR BAM sorting needs many fds) ───────────────────
ulimit -n 65536

# ── Full-genome group (v1_a builds the shared index) ─────────────────────────

echo ""
echo "━━━━ v1_a : full genome | exon-flank | PE | 100 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v1_a.sh"

echo ""
echo "━━━━ v1_b : full genome | exon-flank | SE | 100 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v1_b.sh"

echo ""
echo "━━━━ v1_c : full genome | MANE RNA   | SE | 100 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v1_c.sh"

echo ""
echo "━━━━ v3_a : full genome | exon-flank | PE |  75 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v3_a.sh"

echo ""
echo "━━━━ v4_a : full genome | exon-flank | PE | 150 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v4_a.sh"

# ── Dedup-CDS group (v2_a builds the shared CDS index) ───────────────────────

echo ""
echo "━━━━ v2_a : dedup CDS | exon-flank | PE | 100 bp ━━━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_a.sh"

echo ""
echo "━━━━ v2_b : dedup CDS | exon-flank | SE | 100 bp ━━━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_b.sh"

echo ""
echo "━━━━ v2_c : dedup CDS | MANE RNA   | SE | 100 bp ━━━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_c.sh"

echo ""
echo "━━━━ v2_d : dedup CDS | MANE RNA   | SE | 100 bp | alignIntronMax 1 ━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_d.sh"

# ── MANE transcriptome group (v5_a builds the shared MANE index) ─────────────

echo ""
echo "━━━━ v5_a : MANE ref  | MANE RNA    | SE | 100 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v5_a.sh"

echo ""
echo "━━━━ v5_b : MANE ref  | exon-flank  | PE | 100 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v5_b.sh"

echo ""
echo "━━━━ v5_c : MANE ref  | MANE RNA    | SE |  75 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v5_c.sh"

echo ""
echo "━━━━ v5_d : MANE ref  | MANE RNA    | SE | 200 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v5_d.sh"

# ── Summary ───────────────────────────────────────────────────────────────────

echo ""
echo "════════════════════════════════════════════════════════"
echo "All scripts complete: $(date)"
echo ""
echo "Result TSVs:"
find "$SCRIPT_DIR/results" -name "transcript_uniqueness_factors*.tsv" | sort
echo "════════════════════════════════════════════════════════"
