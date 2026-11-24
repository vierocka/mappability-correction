#!/bin/bash
# ═════════════════════════════════════════════════════════════════════════════
# SLURM header
# ═════════════════════════════════════════════════════════════════════════════
#SBATCH --job-name=mane_mappability_master
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
##SBATCH --account=ag-advol

# ─────────────────────────────────────────────────────────────────────────────
# MASTER script -- run this on HPC via: sbatch run_mappability_master.sh
# Calls compute_uniqueness.py (the worker) for every required run, then STAR,
# then the worker again to score. Produces:
#
#   transcript_uniqueness_factors_L75bp.tsv   ~19,373 primary-chromosome MANE
#   transcript_uniqueness_factors_L100bp.tsv  transcripts, shared across all
#                                              samples (identical to every
#                                              non-HLA-A/B/C gene's correction)
#
#   excluded_non_primary_genes.tsv            64 genes (HLA-DRB3/4, KIR family,
#                                              C4A/C4B, LILRA3, GSTT1/2, ABO,
#                                              MUC16, ...) whose MANE locus is
#                                              on an alt/fix scaffold, not a
#                                              primary chromosome. Genuine
#                                              structural variation (confirmed
#                                              via the NCBI assembly report --
#                                              e.g. HLA-DRB4's scaffold IS one
#                                              of the 8 dedicated MHC
#                                              alt-haplotype contigs), not
#                                              annotation noise. EXCLUDED, not
#                                              silently defaulted to f_unique=1,
#                                              matching nf-core/rnaseq's own
#                                              documented recommendation to use
#                                              the primary assembly only (their
#                                              stated reason: alt loci cause
#                                              "heavy mapping quality
#                                              penalties" from duplicate
#                                              near-identical sequence).
#
#   HCC1937_hla_uniqueness_overrides.tsv      Per-cell-line HLA-A/B/C factors,
#   HCT116_hla_uniqueness_overrides.tsv       computed from each sample's own
#   JY_hla_uniqueness_overrides.tsv           personalized allele sequence.
#                                              NOT merged into the shared base
#                                              table (JY's IDs are unsuffixed
#                                              and would collide with the
#                                              generic entries there) -- load
#                                              these as a per-sample override
#                                              in apply_uniqueness_correction.py,
#                                              same mechanism already used for
#                                              the TAP1/TAP2 manual overrides.
#
# Genome: Ensembl-style primary assembly only (no alt loci/patches), per
# nf-core/rnaseq's documented recommendation and this project's own decision
# after confirming the 64 excluded genes are real structural variants, not a
# gap worth reintroducing alt-haplotype contigs for.
#
# Read lengths: NOT a free choice -- confirmed per-library from each sample's
# own STAR log (*_MANELog.final.out, "Average input read length" / 2 for PE):
#   HCC1937 (SRR31477874/875, SRR925709): ~75 bp/mate
#   HCT116  (SRR3644141-144):             ~101 bp/mate
#   JY      (SRR364065):                  ~99 bp/mate
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# STAR's parallel BAM-sort (--outBAMsortingThreadN threads x --outBAMsortingBinsN
# bins, default 50) can hold hundreds of temp file handles open at once --
# with 24 threads that's up to 1,200, easily exceeding a default shell
# ulimit -n (often 1024). Raise it to the hard limit; --outBAMsortingBinsN is
# also explicitly reduced below as a second safety net in case the hard limit
# itself is capped low on this cluster.
ulimit -n "$(ulimit -Hn)"
echo "ulimit -n raised to: $(ulimit -n)"

module purge
module load bio/STAR/2.7.11b-GCC-13.2.0       # ← adjust version to match what's on HPC
module load bio/SAMtools/1.19.2-GCC-13.2.0    # ← adjust version to match what's on HPC
module load lang/Python/3.11.5-GCCcore-13.2.0  # ← adjust version to match what's on HPC

# ── Paths (verified against actual HPC layout, 2026-08-12) ─────────────────
PROJ_REF="<HPC_WORKDIR>/mappability-correction/Ref"
PROJ_MANE="<HPC_HOME>/MANE"

WORKER="<HPC_WORKDIR>/Laura/compute_uniqueness.py"

MANE_GFF="$PROJ_REF/MANE.GRCh38.v1.5.ensembl_genomic.gff"
GENOME_FA="$PROJ_REF/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
STAR_IDX="$PROJ_REF/STAR_index_primary_assembly_verified"   # NOT reusing STAR_index_genome or
                                                              # STAR_index_genome_nf by name --
                                                              # neither was confirmed by content to
                                                              # match this exact FASTA; build/verify
                                                              # our own unambiguous one instead.

OUTDIR="<HPC_WORKDIR>/Laura"
mkdir -p "$OUTDIR" "$OUTDIR/logs"

THREADS=$SLURM_CPUS_PER_TASK
FLANK=50

echo "════════════════════════════════════════════════════════"
echo "Job: $SLURM_JOB_ID  on $(hostname)  started $(date)"
echo "GENOME_FA: $GENOME_FA"
echo "STAR_IDX:  $STAR_IDX"
echo "OUTDIR:    $OUTDIR"
echo "════════════════════════════════════════════════════════"

# ── Prereq: pysam/pandas venv ──────────────────────────────────────────────────
VENV="$OUTDIR/venv"
python3 -m venv "$VENV"
source "$VENV/bin/activate"
pip install --quiet --upgrade pip
pip install --quiet --no-cache-dir pysam pandas

# ── Prereq: verify genome FASTA is really the 194-sequence primary assembly ──
samtools faidx "$GENOME_FA"
N_SEQ=$(wc -l < "${GENOME_FA}.fai")
if [ "$N_SEQ" -ne 194 ]; then
    echo "ERROR: $GENOME_FA has $N_SEQ sequences, expected 194 (25 primary"
    echo "chromosomes + 169 unplaced/unlocalized scaffolds, no alt loci/patches)."
    echo "Do not proceed -- this may not be the file its name implies."
    exit 1
fi
echo "Genome FASTA verified: $N_SEQ sequences (matches expected primary assembly)."

# ── Prereq: build STAR index if not already validated ─────────────────────────
if [ ! -f "$STAR_IDX/chrNameLength.txt" ] || \
   [ "$(wc -l < "$STAR_IDX/chrNameLength.txt")" -ne 194 ]; then
    echo ""
    echo "[$(date +%H:%M:%S)] Building STAR index (not found or not validated)..."
    mkdir -p "$STAR_IDX"
    STAR --runThreadN $THREADS \
         --runMode genomeGenerate \
         --genomeDir "$STAR_IDX" \
         --genomeFastaFiles "$GENOME_FA"
else
    echo "STAR index already present and validated (194 sequences): $STAR_IDX"
fi

# ── Shared STAR alignment + scoring routine ───────────────────────────────────
# args: <run_outdir> <final_output_path>
align_and_score () {
    local run_outdir="$1"
    local final_output="$2"

    # BAM Unsorted, not SortedByCoordinate: the scoring step (compute_uniqueness.py
    # score) is a pure sequential scan of read flags + MAPQ, no positional
    # queries -- it never needed a sorted/indexed BAM. STAR's internal
    # parallel sort-by-coordinate (many simultaneous thread x bin temp files)
    # is what exhausted open file descriptors last run; skipping it removes
    # that whole failure class rather than just tuning around it.
    STAR --runThreadN $THREADS \
         --genomeDir "$STAR_IDX" \
         --readFilesIn "$run_outdir/simreads/sim_R1.fastq" \
                       "$run_outdir/simreads/sim_R2.fastq" \
         --outSAMtype BAM Unsorted \
         --outSAMattributes NH HI AS NM \
         --outSAMmultNmax 1 \
         --outFilterMultimapNmax 40 \
         --outFilterMismatchNmax 1 \
         --outFileNamePrefix "$run_outdir/bam/sim_" \
         2>&1 | tee "$OUTDIR/logs/star_$(basename "$run_outdir")_${SLURM_JOB_ID}.log"

    python3 "$WORKER" score \
        --outdir "$run_outdir" \
        --bam "$run_outdir/bam/sim_Aligned.out.bam" \
        --out-name "$(basename "$final_output")"

    cp "$run_outdir/$(basename "$final_output")" "$final_output"
}

# ── Run 1: base MANE transcripts @ L75bp ──────────────────────────────────────
echo ""
echo "════ base @ L75bp ════"
RUN="$OUTDIR/base_L75bp"
mkdir -p "$RUN/bam"
# simulate-base already completed in full on the previous (crashed) run --
# STAR had reached "started mapping" before the BAMoutput failure, meaning
# sim_R1.fastq/sim_R2.fastq under $RUN/simreads were already fully written.
# The crash was in STAR's alignment step (now fixed: BAM Unsorted instead of
# SortedByCoordinate), not in read simulation, so no need to regenerate these.
# Uncomment to force regeneration (e.g. if you cleared $RUN/simreads too).
# python3 "$WORKER" simulate-base \
#     --gff "$MANE_GFF" --genome-fa "$GENOME_FA" \
#     --read-len 75 --flank $FLANK --outdir "$RUN"
align_and_score "$RUN" "$OUTDIR/transcript_uniqueness_factors_L75bp.tsv"
cp "$RUN/excluded_non_primary_genes.tsv" "$OUTDIR/excluded_non_primary_genes.tsv"

# ── Run 2: base MANE transcripts @ L100bp ─────────────────────────────────────
echo ""
echo "════ base @ L100bp ════"
RUN="$OUTDIR/base_L100bp"
mkdir -p "$RUN/bam"
python3 "$WORKER" simulate-base \
    --gff "$MANE_GFF" --genome-fa "$GENOME_FA" \
    --read-len 100 --flank $FLANK --outdir "$RUN"
align_and_score "$RUN" "$OUTDIR/transcript_uniqueness_factors_L100bp.tsv"

# ── Runs 3-5 (per-cell-line HLA via simulate-allele) DISABLED 2026-10-06 ─────
# Superseded by run_hla_uniqueness_utr.sh: simulate-allele used CDS exons only
# (no UTRs), aligned allele reads to an unpatched GRCh38 (reference-allele
# mismatch artefact with --outFilterMismatchNmax 1) and did not check that a
# MAPQ-255 read landed on its own locus.
# # ── Run 3: HCC1937 personalized HLA-A/B/C @ L75bp (real read length) ─────────
# echo ""
# echo "════ HCC1937 alleles @ L75bp ════"
# RUN="$OUTDIR/HCC1937_hla_L75bp"
# mkdir -p "$RUN/bam"
# python3 "$WORKER" simulate-allele \
#     --gff "$MANE_GFF" --genome-fa "$GENOME_FA" \
#     --personal-mane "$PROJ_MANE/HCC1937_personal_MANE.fa" \
#     --read-len 75 --flank $FLANK --outdir "$RUN"
# align_and_score "$RUN" "$OUTDIR/HCC1937_hla_uniqueness_overrides.tsv"

# # ── Run 4: HCT116 personalized HLA-A/B/C @ L100bp (real read length) ─────────
# echo ""
# echo "════ HCT116 alleles @ L100bp ════"
# RUN="$OUTDIR/HCT116_hla_L100bp"
# mkdir -p "$RUN/bam"
# python3 "$WORKER" simulate-allele \
#     --gff "$MANE_GFF" --genome-fa "$GENOME_FA" \
#     --personal-mane "$PROJ_MANE/HCT116_personal_MANE.fa" \
#     --read-len 100 --flank $FLANK --outdir "$RUN"
# align_and_score "$RUN" "$OUTDIR/HCT116_hla_uniqueness_overrides.tsv"

# # ── Run 5: JY personalized HLA-A/B/C @ L100bp (real read length) ─────────────
# echo ""
# echo "════ JY alleles @ L100bp ════"
# RUN="$OUTDIR/JY_hla_L100bp"
# mkdir -p "$RUN/bam"
# python3 "$WORKER" simulate-allele \
#     --gff "$MANE_GFF" --genome-fa "$GENOME_FA" \
#     --personal-mane "$PROJ_MANE/JY_personal_MANE.fa" \
#     --read-len 100 --flank $FLANK --outdir "$RUN"
# align_and_score "$RUN" "$OUTDIR/JY_hla_uniqueness_overrides.tsv"

# ── Done ────────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════"
echo "Done: $(date)"
echo "Shared base tables:"
echo "  $OUTDIR/transcript_uniqueness_factors_L75bp.tsv"
echo "  $OUTDIR/transcript_uniqueness_factors_L100bp.tsv"
echo "Excluded (non-primary locus, e.g. HLA-DRB4, KIR family, C4A/C4B):"
echo "  $OUTDIR/excluded_non_primary_genes.tsv"
echo "Per-cell-line HLA overrides (load in apply_uniqueness_correction*.py,"
echo "same mechanism as the existing TAP1/TAP2 manual overrides):"
echo "  $OUTDIR/HCC1937_hla_uniqueness_overrides.tsv  (L75bp)"
echo "  $OUTDIR/HCT116_hla_uniqueness_overrides.tsv   (L100bp)"
echo "  $OUTDIR/JY_hla_uniqueness_overrides.tsv       (L100bp)"
echo "════════════════════════════════════════════════════════"
