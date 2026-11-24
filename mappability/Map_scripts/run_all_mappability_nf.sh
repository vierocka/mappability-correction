#!/bin/bash
#SBATCH --cpus-per-task=64
#SBATCH --mem=192gb
#SBATCH --time=72:00:00
#SBATCH --output=mappability_all_nf_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --account=your_slurm_account

# ─────────────────────────────────────────────────────────────────────────────
# nf-reference (Ensembl GRCh38.116) mappability — HPC submission script
#
# Runs all 42 "_nf" scripts: the genomic-reference arm (v1_a-d, v3_a-d,
# v4_a-d, v6_a-b, v8_a — 15 scripts), the genome-wide dedup-CDS arm (v7_a-j,
# v8_b, v8_c — 12 scripts), and the MANE-restricted dedup-CDS arm (v2_a-o —
# 15 scripts, a third, separate CDS reference from the genome-wide one), in
# an order that guarantees every shared build step happens exactly once
# before anything that depends on it:
#   STAR_index_genome_nf        — built by v1_a_nf
#   MANE_RNA_nf.fna             — built by v1_c_nf
#   dedup_cds_nf.fna            — built by v7_a_nf (genome-wide protein-coding
#                                  CDS extracted from the nf GTF, then seqkit
#                                  rmdup'd — see mappability_correction_v7_a_nf.sh)
#   STAR_index_dedup_cds_nf     — built by v7_a_nf
#   STAR_index_MANE_nf          — built by v8_c_nf
#   dedup_cds_mane_nf.fna       — built by v2_a_nf (MANE_Select-only CDS,
#                                  same splicing approach, smaller scope —
#                                  see mappability_correction_v2_a_nf.sh)
#   STAR_index_dedup_cds_mane_nf — built by v2_a_nf
# Every other script reuses what's already built (see each script's own
# Step 0/0b/0c/0d/0e skip-if-exists guards). v8_b_nf/v8_c_nf additionally
# require v8_a_nf's simulated reads to already exist (checked explicitly;
# they abort with a clear error otherwise).
#
# seqkit (needed once each, by v7_a_nf's and v2_a_nf's Step 0d) has no
# confirmed module name on this cluster — if either step fails, load seqkit
# manually first and re-run, or dedup locally and scp the result back (see
# each script's Step 0d instructions).
#
# Skips any script whose final TSV already exists, so a partially completed
# run can be resubmitted safely without redoing finished work.
#
# Does not touch the existing (non-_nf) pipeline, its scripts, indices, or
# results in any way.
#
# Output: one TSV per script in results/mappability_*_nf/
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# Resolved from the script's own location (Map_scripts/ -> up two levels to
# mappability-correction/) rather than a hardcoded absolute path, so this
# works regardless of where the repo is cloned/deployed.
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
MAP="$BASE/mappability/Map_scripts"
RES="$BASE/mappability/results"

mkdir -p "$RES"

module purge
module load bio/STAR/2.7.11b-GCC-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0
module load lang/Python/3.11.5-GCCcore-13.2.0

ulimit -n 65536

echo "════════════════════════════════════════════════════════"
echo "nf-reference mappability — full run (42 scripts)"
echo "Started : $(date)"
echo "Host    : $(hostname)"
echo "Base    : $BASE"
echo "════════════════════════════════════════════════════════"

run() {
    local SID="$1"
    local OUTDIR_NAME="$2"
    local TSV_NAME="$3"

    local OUT_TSV="$RES/$OUTDIR_NAME/$TSV_NAME"

    if [ -f "$OUT_TSV" ]; then
        echo "SKIP $SID — already done: $OUT_TSV"
        return
    fi

    echo ""
    echo "━━━━ $SID ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    bash "$MAP/mappability_correction_${SID}.sh"
}

# ── Prerequisite pair — builds the two shared resources (STAR index, MANE_RNA_nf.fna) ──
run v1_a_nf mappability_genomic_nf                       transcript_uniqueness_factors_genomic_L100bp_nf.tsv
run v1_c_nf mappability_genomic_RNA_SE_L100bp_nf         transcript_uniqueness_factors_genomic_RNA_SE_L100bp_nf.tsv

# ── v1 — full genome, 100 bp (remaining) ─────────────────────────────────────
run v1_b_nf mappability_genomic_SE_L100bp_nf             transcript_uniqueness_factors_genomic_SE_L100bp_nf.tsv
run v1_d_nf mappability_genomic_RNA_PE_nf                transcript_uniqueness_factors_genomic_RNA_PE_nf.tsv

# ── v3 — full genome, 75 bp ───────────────────────────────────────────────────
run v3_a_nf mappability_genomic_L75bp_nf                 transcript_uniqueness_factors_genomic_L75bp_nf.tsv
run v3_b_nf mappability_genomic_SE_L75bp_nf              transcript_uniqueness_factors_genomic_SE_L75bp_nf.tsv
run v3_c_nf mappability_genomic_RNA_SE_L75bp_nf          transcript_uniqueness_factors_genomic_RNA_SE_L75bp_nf.tsv
run v3_d_nf mappability_genomic_RNA_PE_L75bp_nf          transcript_uniqueness_factors_genomic_RNA_PE_L75bp_nf.tsv

# ── v4 — full genome, 150 bp ──────────────────────────────────────────────────
run v4_a_nf mappability_genomic_L150bp_nf                transcript_uniqueness_factors_genomic_L150bp_nf.tsv
run v4_b_nf mappability_genomic_SE_L150bp_nf             transcript_uniqueness_factors_genomic_SE_L150bp_nf.tsv
run v4_c_nf mappability_genomic_RNA_SE_L150bp_nf         transcript_uniqueness_factors_genomic_RNA_SE_L150bp_nf.tsv
run v4_d_nf mappability_genomic_RNA_PE_L150bp_nf         transcript_uniqueness_factors_genomic_RNA_PE_L150bp_nf.tsv

# ── v6 — full genome, 200 bp ──────────────────────────────────────────────────
run v6_a_nf mappability_genomic_RNA_SE_L200bp_nf         transcript_uniqueness_factors_genomic_RNA_SE_L200bp_nf.tsv
run v6_b_nf mappability_genomic_RNA_PE_L200bp_nf         transcript_uniqueness_factors_genomic_RNA_PE_L200bp_nf.tsv

# ── v8 — proper PE 150 bp (realistic library geometry) ───────────────────────
run v8_a_nf mappability_genomic_RNA_PE_L150bp_proper_nf  transcript_uniqueness_factors_genomic_RNA_PE_L150bp_proper_nf.tsv

# ── v7 — dedup-CDS reads (genome-wide, protein-coding), full genome target ───
# v7_a_nf builds dedup_cds_nf.fna + STAR_index_dedup_cds_nf (shared, once)
run v7_a_nf mappability_dedup_cds_reads_SE_nf            cds_uniqueness_factors_genome_SE_nf.tsv
run v7_b_nf mappability_dedup_cds_reads_PE_nf            cds_uniqueness_factors_genome_PE_nf.tsv
run v7_c_nf mappability_dedup_cds_reads_SE_L75bp_nf      cds_uniqueness_factors_genome_SE_L75bp_nf.tsv
run v7_d_nf mappability_dedup_cds_reads_PE_L75bp_nf      cds_uniqueness_factors_genome_PE_L75bp_nf.tsv
run v7_e_nf mappability_dedup_cds_reads_SE_L150bp_nf     cds_uniqueness_factors_genome_SE_L150bp_nf.tsv
run v7_f_nf mappability_dedup_cds_reads_PE_L150bp_nf     cds_uniqueness_factors_genome_PE_L150bp_nf.tsv
run v7_i_nf mappability_dedup_cds_reads_SE_L200bp_nf     cds_uniqueness_factors_genome_SE_L200bp_nf.tsv
run v7_j_nf mappability_dedup_cds_reads_PE_L200bp_nf     cds_uniqueness_factors_genome_PE_L200bp_nf.tsv

# ── v7 — dedup-CDS reads, dedup-CDS self-map (residual within-CDS redundancy) ─
run v7_g_nf mappability_dedup_cds_reads_selfmap_SE_nf    cds_uniqueness_factors_selfmap_SE_nf.tsv
run v7_h_nf mappability_dedup_cds_reads_selfmap_PE_nf    cds_uniqueness_factors_selfmap_PE_nf.tsv

# ── v8_b/c — MANE RNA proper-PE reads (reused from v8_a_nf), dedup-CDS / MANE-self targets ──
run v8_b_nf mappability_dedup_cds_RNA_PE_L150bp_proper_nf     transcript_uniqueness_factors_dedup_cds_RNA_PE_L150bp_proper_nf.tsv
run v8_c_nf mappability_MANE_RNA_PE_L150bp_proper_nf          transcript_uniqueness_factors_MANE_RNA_PE_L150bp_proper_nf.tsv

# ── v2 — MANE-restricted dedup-CDS reference (separate from v7's genome-wide one) ──
# v2_a_nf builds dedup_cds_mane_nf.fna + STAR_index_dedup_cds_mane_nf (shared, once)
run v2_a_nf mappability_dedup_cds_PE_flank_nf             transcript_uniqueness_factors_dedup_cds_PE_flank_L100bp_nf.tsv
run v2_b_nf mappability_dedup_cds_SE_flank_nf             transcript_uniqueness_factors_dedup_cds_SE_flank_L100bp_nf.tsv
run v2_c_nf mappability_dedup_cds_nf                      transcript_uniqueness_factors_dedup_cds_L100bp_nf.tsv
run v2_d_nf mappability_dedup_cds_alignIntronMax1_nf      transcript_uniqueness_factors_dedup_cds_alignIntronMax1_L100bp_nf.tsv
run v2_e_nf mappability_dedup_cds_RNA_PE_nf               transcript_uniqueness_factors_dedup_cds_RNA_PE_nf.tsv
run v2_f_nf mappability_dedup_cds_PE_flank_L75bp_nf       transcript_uniqueness_factors_dedup_cds_PE_flank_L75bp_nf.tsv
run v2_g_nf mappability_dedup_cds_SE_flank_L75bp_nf       transcript_uniqueness_factors_dedup_cds_SE_flank_L75bp_nf.tsv
run v2_h_nf mappability_dedup_cds_PE_flank_L150bp_nf      transcript_uniqueness_factors_dedup_cds_PE_flank_L150bp_nf.tsv
run v2_i_nf mappability_dedup_cds_SE_flank_L150bp_nf      transcript_uniqueness_factors_dedup_cds_SE_flank_L150bp_nf.tsv
run v2_j_nf mappability_dedup_cds_RNA_SE_L75bp_nf         transcript_uniqueness_factors_dedup_cds_RNA_SE_L75bp_nf.tsv
run v2_k_nf mappability_dedup_cds_RNA_PE_L75bp_nf         transcript_uniqueness_factors_dedup_cds_RNA_PE_L75bp_nf.tsv
run v2_l_nf mappability_dedup_cds_RNA_SE_L150bp_nf        transcript_uniqueness_factors_dedup_cds_RNA_SE_L150bp_nf.tsv
run v2_m_nf mappability_dedup_cds_RNA_PE_L150bp_nf        transcript_uniqueness_factors_dedup_cds_RNA_PE_L150bp_nf.tsv
run v2_n_nf mappability_dedup_cds_RNA_SE_L200bp_nf        transcript_uniqueness_factors_dedup_cds_RNA_SE_L200bp_nf.tsv
run v2_o_nf mappability_dedup_cds_RNA_PE_L200bp_nf        transcript_uniqueness_factors_dedup_cds_RNA_PE_L200bp_nf.tsv

echo ""
echo "════════════════════════════════════════════════════════"
echo "All done: $(date)"
echo "Output TSVs:"
ls -lh "$RES"/mappability_*_nf/transcript_uniqueness_factors_*.tsv 2>/dev/null | awk '{print $5, $9}' | sort
echo "════════════════════════════════════════════════════════"
