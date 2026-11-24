#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH --output=sc3prime_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL

# ─────────────────────────────────────────────────────────────────────────────
# 3' window uniqueness factor analysis — HPC submission script
#
# Computes per-transcript UF restricted to the last 250/400/500/600 bp of
# each transcript for all MANE RNA and exon-flank read settings (v1-v8).
#
# Motivation: 10x Genomics 3' scRNA-seq captures only the last ~200-600 bp
# of each transcript. Genes may have good overall UF but poor 3' end UF,
# making them invisible or noisy in 10x data. This analysis identifies
# which cell marker genes are most affected.
#
# Output: one TSV per run in sc3UTRs/results/
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

BASE="./mappability-correction"
REF="$BASE/Ref"
RES="$BASE/mappability/results"
SC3="$BASE/mappability/sc3UTRs"
OUT="$SC3/results"
PY_UF="$SC3/compute_3prime_uf.py"
PY_SIM="$SC3/compute_3end_similarity.py"

# ── Reference annotations ─────────────────────────────────────────────────────
# GENCODE v49: transcript FASTA already available, no extraction needed.
# GRCh38.p14 assembly — same underlying genome as the STAR index.
# Note: gffread cannot parse the NCBI RefSeq GTF (GRCh38p14_primary.gtf)
# due to non-standard attribute formatting. RefSeq mode is not used.
#
GENCODE_FA="$REF/gencode.v49.transcripts.fa.gz"
GENCODE_GTF="$REF/gencode.v49.basic.annotation.gtf"

mkdir -p "$OUT"

module purge
module load lang/Python/3.11.5-GCCcore-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0

echo "════════════════════════════════════════════════════════"
echo "3' window analyses — MANE RNA BAMs + 3' end similarity"
echo "Started : $(date)"
echo "Host    : $(hostname)"
echo "════════════════════════════════════════════════════════"

# ═══════════════════════════════════════════════════════════════
# PART A: per-BAM 3' window UF (compute_3prime_uf.py)
# ═══════════════════════════════════════════════════════════════

run() {
    local BAM_DIR="$1"
    local TSV_NAME="$2"
    local READ_LEN="$3"
    local LABEL="$4"

    local BAM="$RES/$BAM_DIR/bam/sim_Aligned.sortedByCoord.out.bam"
    local TSV="$RES/$TSV_NAME"
    local OUT_TSV="$OUT/sc3prime_${LABEL}.tsv"

    if [ ! -f "$BAM" ]; then
        echo "SKIP $LABEL — BAM not found: $BAM"
        return
    fi
    if [ ! -f "$TSV" ]; then
        echo "SKIP $LABEL — ref TSV not found: $TSV"
        return
    fi
    if [ -f "$OUT_TSV" ]; then
        echo "SKIP $LABEL — already done: $OUT_TSV"
        return
    fi

    echo ""
    echo "━━━━ $LABEL ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    python3 "$PY_UF" "$BAM" "$TSV" "$READ_LEN" "$OUT_TSV" "$LABEL"
}

# ── v1 — full genome, 100 bp ──────────────────────────────────────────────────
run mappability_genomic_SE           transcript_uniqueness_factors_genomic_SE_L100bp.tsv     100 v1b_genome_exflank_SE_100
run mappability_genomic              transcript_uniqueness_factors_genomic_L100bp.tsv         100 v1a_genome_exflank_PE_100
run mappability_genomic_RNA_SE       transcript_uniqueness_factors_genomic_RNA_SE_L100bp.tsv  100 v1c_genome_RNA_SE_100
run mappability_genomic_RNA_PE       transcript_uniqueness_factors_genomic_RNA_PE.tsv         100 v1d_genome_RNA_PE_100

# ── v3 — full genome, 75 bp (closest to 10x v3 R2 read length ~90 bp) ────────
run mappability_genomic_SE_L75bp     transcript_uniqueness_factors_genomic_SE_L75bp.tsv       75  v3b_genome_exflank_SE_75
run mappability_genomic_L75bp        transcript_uniqueness_factors_genomic_L75bp.tsv           75  v3a_genome_exflank_PE_75
run mappability_genomic_RNA_SE_L75bp transcript_uniqueness_factors_genomic_RNA_SE_L75bp.tsv   75  v3c_genome_RNA_SE_75
run mappability_genomic_RNA_PE_L75bp transcript_uniqueness_factors_genomic_RNA_PE_L75bp.tsv   75  v3d_genome_RNA_PE_75

# ── v4 — full genome, 150 bp ─────────────────────────────────────────────────
run mappability_genomic_SE_L150bp     transcript_uniqueness_factors_genomic_SE_L150bp.tsv     150 v4b_genome_exflank_SE_150
run mappability_genomic_L150bp        transcript_uniqueness_factors_genomic_L150bp.tsv         150 v4a_genome_exflank_PE_150
run mappability_genomic_RNA_SE_L150bp transcript_uniqueness_factors_genomic_RNA_SE_L150bp.tsv 150 v4c_genome_RNA_SE_150
run mappability_genomic_RNA_PE_L150bp transcript_uniqueness_factors_genomic_RNA_PE_L150bp.tsv 150 v4d_genome_RNA_PE_150

# ── v5 — MANE transcriptome reference (upper bound / self-map) ────────────────
run mappability_MANE_RNA_SE_L75bp    transcript_uniqueness_factors_MANE_RNA_SE_L75bp.tsv      75  v5c_MANE_RNA_SE_75
run mappability_MANE_RNA_SE          transcript_uniqueness_factors_MANE_RNA_SE_L100bp.tsv     100 v5a_MANE_RNA_SE_100
run mappability_MANE_RNA_SE_L150bp   transcript_uniqueness_factors_MANE_RNA_SE_L150bp.tsv    150 v5e_MANE_RNA_SE_150
run mappability_MANE_RNA_SE_L200bp   transcript_uniqueness_factors_MANE_RNA_SE_L200bp.tsv    200 v5d_MANE_RNA_SE_200
run mappability_MANE_RNA_PE_L75bp    transcript_uniqueness_factors_MANE_RNA_PE_L75bp.tsv      75  v5f_MANE_RNA_PE_75
run mappability_MANE_RNA_PE          transcript_uniqueness_factors_MANE_RNA_PE.tsv            100 v5g_MANE_RNA_PE_100
run mappability_MANE_RNA_PE_L150bp   transcript_uniqueness_factors_MANE_RNA_PE_L150bp.tsv    150 v5h_MANE_RNA_PE_150
run mappability_MANE_RNA_PE_L200bp   transcript_uniqueness_factors_MANE_RNA_PE_L200bp.tsv    200 v5i_MANE_RNA_PE_200
run mappability_MANE_exonflank_SE_L75bp   transcript_uniqueness_factors_MANE_exonflank_SE_L75bp.tsv  75  v5j_MANE_exflank_SE_75
run mappability_MANE_exonflank_SE         transcript_uniqueness_factors_MANE_exonflank_SE.tsv        100 v5k_MANE_exflank_SE_100
run mappability_MANE_exonflank_SE_L150bp  transcript_uniqueness_factors_MANE_exonflank_SE_L150bp.tsv 150 v5l_MANE_exflank_SE_150
run mappability_MANE_exonflank_PE_L75bp   transcript_uniqueness_factors_MANE_exonflank_PE_L75bp.tsv  75  v5m_MANE_exflank_PE_75
run mappability_MANE_exonflank_PE         transcript_uniqueness_factors_MANE_exonflank_PE_L100bp.tsv 100 v5b_MANE_exflank_PE_100
run mappability_MANE_exonflank_PE_L150bp  transcript_uniqueness_factors_MANE_exonflank_PE_L150bp.tsv 150 v5n_MANE_exflank_PE_150

# ── v6 — full genome, 200 bp ─────────────────────────────────────────────────
run mappability_genomic_RNA_SE_L200bp transcript_uniqueness_factors_genomic_RNA_SE_L200bp.tsv 200 v6a_genome_RNA_SE_200
run mappability_genomic_RNA_PE_L200bp transcript_uniqueness_factors_genomic_RNA_PE_L200bp.tsv 200 v6b_genome_RNA_PE_200

# ── v8 — proper PE 150 bp (most realistic library geometry) ──────────────────
run mappability_genomic_RNA_PE_L150bp_proper     transcript_uniqueness_factors_genomic_RNA_PE_L150bp_proper.tsv  150 v8a_genome_RNA_PE_150_proper
run mappability_dedup_cds_RNA_PE_L150bp_proper   transcript_uniqueness_factors_dedup_cds_RNA_PE_L150bp_proper.tsv 150 v8b_dedupcds_RNA_PE_150_proper
run mappability_MANE_RNA_PE_L150bp_proper        transcript_uniqueness_factors_MANE_RNA_PE_L150bp_proper.tsv     150 v8c_MANE_RNA_PE_150_proper

# ── v2 — dedup CDS reference (protein-coding scope) ──────────────────────────
run mappability_dedup_cds            transcript_uniqueness_factors_dedup_cds_L100bp.tsv         100 v2c_dedupcds_RNA_SE_100
run mappability_dedup_cds_RNA_PE     transcript_uniqueness_factors_dedup_cds_RNA_PE.tsv          100 v2e_dedupcds_RNA_PE_100
run mappability_dedup_cds_RNA_SE_L75bp  transcript_uniqueness_factors_dedup_cds_RNA_SE_L75bp.tsv  75  v2j_dedupcds_RNA_SE_75
run mappability_dedup_cds_RNA_PE_L75bp  transcript_uniqueness_factors_dedup_cds_RNA_PE_L75bp.tsv  75  v2k_dedupcds_RNA_PE_75
run mappability_dedup_cds_RNA_SE_L150bp transcript_uniqueness_factors_dedup_cds_RNA_SE_L150bp.tsv 150 v2l_dedupcds_RNA_SE_150
run mappability_dedup_cds_RNA_PE_L150bp transcript_uniqueness_factors_dedup_cds_RNA_PE_L150bp.tsv 150 v2m_dedupcds_RNA_PE_150
run mappability_dedup_cds_RNA_SE_L200bp transcript_uniqueness_factors_dedup_cds_RNA_SE_L200bp.tsv 200 v2n_dedupcds_RNA_SE_200
run mappability_dedup_cds_RNA_PE_L200bp transcript_uniqueness_factors_dedup_cds_RNA_PE_L200bp.tsv 200 v2o_dedupcds_RNA_PE_200

echo ""
echo "════════════════════════════════════════════════════════"
echo "Part A complete: $(date)"
echo "Output TSVs:"
ls -lh "$OUT"/sc3prime_*.tsv 2>/dev/null | awk '{print $5, $9}' | sort

# ═══════════════════════════════════════════════════════════════
# PART B: 3' end sequence similarity (compute_3end_similarity.py)
# ═══════════════════════════════════════════════════════════════

echo ""
echo "════════════════════════════════════════════════════════"
echo "Part B: 3' end sequence similarity — GENCODE v49"
echo "════════════════════════════════════════════════════════"

python3 "$PY_SIM" \
    --fasta  "$GENCODE_FA" \
    --gtf    "$GENCODE_GTF" \
    --format gencode \
    --label  gencode_v49 \
    --outdir "$OUT"

echo ""
echo "════════════════════════════════════════════════════════"
echo "All done: $(date)"
echo "════════════════════════════════════════════════════════"
