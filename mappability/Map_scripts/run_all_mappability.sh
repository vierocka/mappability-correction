#!/bin/bash -l
#SBATCH --cpus-per-task=64
#SBATCH --mem=192gb
#SBATCH --time=72:00:00
#SBATCH --output=mappability_all_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL

# ─────────────────────────────────────────────────────────────────────────────
# Run all mappability correction scripts sequentially (43 STAR + 10 dedup CDS reads = 53 total).
#
# Execution order matters:
#   v1_a runs first  → builds the shared genome STAR index (STAR_index_genome)
#   v2_a runs first of the CDS group → builds STAR_index_dedup_cds
#   All subsequent scripts reuse whichever index already exists.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

echo "════════════════════════════════════════════════════════"
echo "Mappability correction — full run (43 scripts)"
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

# ── DONE — completed successfully; TSV output verified ───────────────────────
# echo ""
# echo "━━━━ v1_a : full genome | exon-flank | PE | 100 bp ━━━━━━━━━━━━━━━━━━━━"
# bash "$SCRIPT_DIR/mappability_correction_v1_a.sh"
#
# echo ""
# echo "━━━━ v1_b : full genome | exon-flank | SE | 100 bp ━━━━━━━━━━━━━━━━━━━━"
# bash "$SCRIPT_DIR/mappability_correction_v1_b.sh"
#
# echo ""
# echo "━━━━ v1_c : full genome | MANE RNA   | SE | 100 bp ━━━━━━━━━━━━━━━━━━━━"
# bash "$SCRIPT_DIR/mappability_correction_v1_c.sh"
#
# echo ""
# echo "━━━━ v3_a : full genome | exon-flank | PE |  75 bp ━━━━━━━━━━━━━━━━━━━━"
# bash "$SCRIPT_DIR/mappability_correction_v3_a.sh"
#
# echo ""
# echo "━━━━ v4_a : full genome | exon-flank | PE | 150 bp ━━━━━━━━━━━━━━━━━━━━"
# bash "$SCRIPT_DIR/mappability_correction_v4_a.sh"

echo ""
echo "━━━━ v1_d : full genome | MANE RNA   | PE | 100 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v1_d.sh"

echo ""
echo "━━━━ v3_b : full genome | exon-flank | SE |  75 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v3_b.sh"

echo ""
echo "━━━━ v3_c : full genome | MANE RNA   | SE |  75 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v3_c.sh"

echo ""
echo "━━━━ v3_d : full genome | MANE RNA   | PE |  75 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v3_d.sh"

echo ""
echo "━━━━ v4_b : full genome | exon-flank | SE | 150 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v4_b.sh"

echo ""
echo "━━━━ v4_c : full genome | MANE RNA   | SE | 150 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v4_c.sh"

echo ""
echo "━━━━ v4_d : full genome | MANE RNA   | PE | 150 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v4_d.sh"

echo ""
echo "━━━━ v6_a : full genome | MANE RNA   | SE | 200 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v6_a.sh"

echo ""
echo "━━━━ v6_b : full genome | MANE RNA   | PE | 200 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v6_b.sh"

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

echo ""
echo "━━━━ v2_e : dedup CDS | MANE RNA   | PE | 100 bp ━━━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_e.sh"

echo ""
echo "━━━━ v2_f : dedup CDS | exon-flank | PE |  75 bp ━━━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_f.sh"

echo ""
echo "━━━━ v2_g : dedup CDS | exon-flank | SE |  75 bp ━━━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_g.sh"

echo ""
echo "━━━━ v2_h : dedup CDS | exon-flank | PE | 150 bp ━━━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_h.sh"

echo ""
echo "━━━━ v2_i : dedup CDS | exon-flank | SE | 150 bp ━━━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_i.sh"

echo ""
echo "━━━━ v2_j : dedup CDS | MANE RNA   | SE |  75 bp ━━━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_j.sh"

echo ""
echo "━━━━ v2_k : dedup CDS | MANE RNA   | PE |  75 bp ━━━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_k.sh"

echo ""
echo "━━━━ v2_l : dedup CDS | MANE RNA   | SE | 150 bp ━━━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_l.sh"

echo ""
echo "━━━━ v2_m : dedup CDS | MANE RNA   | PE | 150 bp ━━━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_m.sh"

echo ""
echo "━━━━ v2_n : dedup CDS | MANE RNA   | SE | 200 bp ━━━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_n.sh"

echo ""
echo "━━━━ v2_o : dedup CDS | MANE RNA   | PE | 200 bp ━━━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v2_o.sh"

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

echo ""
echo "━━━━ v5_e : MANE ref  | MANE RNA    | SE | 150 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v5_e.sh"

echo ""
echo "━━━━ v5_f : MANE ref  | MANE RNA    | PE |  75 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v5_f.sh"

echo ""
echo "━━━━ v5_g : MANE ref  | MANE RNA    | PE | 100 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v5_g.sh"

echo ""
echo "━━━━ v5_h : MANE ref  | MANE RNA    | PE | 150 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v5_h.sh"

echo ""
echo "━━━━ v5_i : MANE ref  | MANE RNA    | PE | 200 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v5_i.sh"

echo ""
echo "━━━━ v5_j : MANE ref  | exon-flank  | SE |  75 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v5_j.sh"

echo ""
echo "━━━━ v5_k : MANE ref  | exon-flank  | SE | 100 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v5_k.sh"

echo ""
echo "━━━━ v5_l : MANE ref  | exon-flank  | SE | 150 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v5_l.sh"

echo ""
echo "━━━━ v5_m : MANE ref  | exon-flank  | PE |  75 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v5_m.sh"

echo ""
echo "━━━━ v5_n : MANE ref  | exon-flank  | PE | 150 bp ━━━━━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v5_n.sh"

# ── Dedup CDS reads group (v7) — reads from CDS FASTA; STAR alignment ─────────
# Run after v1 and v2 groups (genome and CDS STAR indices must exist first).
# v7_a-f,i,j: reads vs full genome; v7_g-h: self-map vs dedup CDS.
# Output directories align with kallisto_mappability_crossval.sh expectations.

echo ""
echo "━━━━ v7_a : dedup CDS reads | full genome | SE | 100 bp ━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v7_a.sh"

echo ""
echo "━━━━ v7_b : dedup CDS reads | full genome | PE | 100 bp ━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v7_b.sh"

echo ""
echo "━━━━ v7_c : dedup CDS reads | full genome | SE |  75 bp ━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v7_c.sh"

echo ""
echo "━━━━ v7_d : dedup CDS reads | full genome | PE |  75 bp ━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v7_d.sh"

echo ""
echo "━━━━ v7_e : dedup CDS reads | full genome | SE | 150 bp ━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v7_e.sh"

echo ""
echo "━━━━ v7_f : dedup CDS reads | full genome | PE | 150 bp ━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v7_f.sh"

echo ""
echo "━━━━ v7_g : dedup CDS reads | dedup CDS   | SE | 100 bp (self-map) ━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v7_g.sh"

echo ""
echo "━━━━ v7_h : dedup CDS reads | dedup CDS   | PE | 100 bp (self-map) ━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v7_h.sh"

echo ""
echo "━━━━ v7_i : dedup CDS reads | full genome | SE | 200 bp ━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v7_i.sh"

echo ""
echo "━━━━ v7_j : dedup CDS reads | full genome | PE | 200 bp ━━━━━━━━━━━━━━━━"
bash "$SCRIPT_DIR/mappability_correction_v7_j.sh"

# ── Summary ───────────────────────────────────────────────────────────────────

echo ""
echo "════════════════════════════════════════════════════════"
echo "All scripts complete: $(date)"
echo ""
echo "Result TSVs:"
find "$SCRIPT_DIR/results" -name "transcript_uniqueness_factors*.tsv" | sort
echo "════════════════════════════════════════════════════════"
