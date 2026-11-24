#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --mem=192gb
#SBATCH --time=24:00:00
#SBATCH --output=MANE_personalized_UTR_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>
#SBATCH --account=<group>

###############################################################################
# MANE-personalized STAR + featureCounts, rerun with UTR-fixed references
# (2026-10-05). Replaces Laura_MANE_personalized.sh.
#
# Changes vs Laura_MANE_personalized.sh:
#   1. References: {JY,HCC1937,HCT116}_personal_MANE.fa / _cds_mane.saf from
#      build_personal_MANE_utr.py -- HLA-A/B/C = MANE 5'UTR + allele CDS +
#      MANE 3'UTR. The old ones were CDS-only for HLA in ALL three samples
#      (JY included), while every other gene is counted over its whole
#      transcript (SAF End = transcript length).
#   2. STAR settings identical for every sample. Old script had
#      --alignIntronMax 1 --alignIntronMin 2 for HCT116 only; JY and HCC1937
#      were aligned with splicing allowed against a transcript reference.
#   3. nf-core/rnaseq-like STAR options where they apply to a transcriptome
#      reference: --outFilterMultimapNmax 20, --outSAMattributes NH HI AS NM MD,
#      --runRNGseed 0. (No GTF / two-pass / TranscriptomeSAM -- not applicable,
#      the "genome" here IS the transcriptome.)
#   4. featureCounts run twice per sample:
#        *_fc.counts.csv   unique reads only (same counting as before)
#        *_fcM.counts.csv  -M --fraction: multimappers split 1/NH
#      Heterozygous HLA alleles (HCC1937, HCT116) are near-identical and now
#      share identical MANE UTRs, so most of their reads are NH=2 and are
#      DROPPED by unique-only counting, whereas JY (homozygous, one entry per
#      gene) keeps them. Use the fcM table and sum the two allele rows per
#      gene for any HLA comparison across cell lines.
#   5. Pre-flight: every reference file and every R1/R2 FASTQ is checked
#      (exists, non-empty) BEFORE anything runs; job aborts listing all
#      missing files. Old STAR indexes are removed and rebuilt (old script's
#      plain mkdir did not clear them).
#   6. New output prefix *_MANEpersUTR so the previous *_MANEpers results
#      stay in place for comparison.
#
# Reads: raw FASTQs, untrimmed (as before). Trimmomatic output exists only for
# HCC1937 (<HPC_WORKDIR>/Laura/hcc1937_trimmomatic) -- trimming must be
# applied to all samples or none to keep them comparable.
###############################################################################

set -euo pipefail

THREADS=24
REF_DIR="<HPC_HOME>/MANE"           # personal_MANE_utr/* copied here
WORK_DIR="<HPC_WORKDIR>/Laura"
TAG="MANEpersUTR"

# cell_line  run_accession  sample_dir (under WORK_DIR)
SAMPLES="JY      SRR364065    SRR364065_JY
HCT116  SRR3644141   SRR3644141_HCT116_total
HCT116  SRR3644142   SRR3644142_HCT116
HCT116  SRR3644143   SRR3644143_HCT116
HCT116  SRR3644144   SRR3644144_HCT116
HCC1937 SRR925709    SRR925709_HCC1937
HCC1937 SRR31477874  SRR31477874_HCC1937
HCC1937 SRR31477875  SRR31477875_HCC1937"

CELL_LINES="JY HCT116 HCC1937"

################# pre-flight: references + FASTQs ################
missing=()
for cl in $CELL_LINES; do
  for f in "${REF_DIR}/${cl}_personal_MANE.fa" "${REF_DIR}/${cl}_cds_mane.saf"; do
    [ -s "$f" ] || missing+=("$f")
  done
  # must be the UTR-fixed reference: HLA headers carry a cds:<start>-<end> tag
  if [ -s "${REF_DIR}/${cl}_personal_MANE.fa" ] && \
     ! grep -q '^>ENST00000376809.* cds:[0-9]*-[0-9]*' "${REF_DIR}/${cl}_personal_MANE.fa"; then
    missing+=("${REF_DIR}/${cl}_personal_MANE.fa is the OLD CDS-only reference (no cds: tag on HLA-A)")
  fi
done

while read -r cl run dir; do
  [ -z "$cl" ] && continue
  for r in 1 2; do
    fq="${WORK_DIR}/${dir}/${run}_${r}.fastq"
    [ -s "$fq" ] || missing+=("$fq")
  done
done <<< "$SAMPLES"

if [ ${#missing[@]} -gt 0 ]; then
  echo "ERROR: missing or empty input(s) -- nothing was run:" >&2
  printf '  %s\n' "${missing[@]}" >&2
  exit 1
fi
echo "Pre-flight OK: all references and FASTQs present."

################# STAR ################
set +u   # module scripts reference unset variables
module unload compiler/GCC 2>/dev/null || true
module add bio/STAR/2.7.10b-GCC-11.3.0
set -u

for cl in $CELL_LINES; do
  idx="${REF_DIR}/STAR_index_MANE_${cl}"
  echo "=== STAR index ${cl} ==="
  rm -rf "$idx"
  mkdir -p "$idx"
  STAR --runThreadN $THREADS \
       --runMode genomeGenerate \
       --genomeDir "$idx" \
       --genomeFastaFiles "${REF_DIR}/${cl}_personal_MANE.fa" \
       --genomeSAindexNbases 12 \
       --genomeChrBinNbits 11 \
       --outFileNamePrefix "${idx}/"
done

while read -r cl run dir; do
  [ -z "$cl" ] && continue
  sdir="${WORK_DIR}/${dir}"
  echo "=== STAR align ${run} (${cl}) ==="
  STAR \
    --runMode alignReads \
    --runThreadN $THREADS \
    --genomeDir "${REF_DIR}/STAR_index_MANE_${cl}" \
    --readFilesIn "${sdir}/${run}_1.fastq" "${sdir}/${run}_2.fastq" \
    --alignIntronMax 1 \
    --alignIntronMin 2 \
    --outFilterMultimapNmax 20 \
    --outSAMattributes NH HI AS NM MD \
    --runRNGseed 0 \
    --outBAMsortingBinsN 20 \
    --limitBAMsortRAM 10000000000 \
    --outFileNamePrefix "${sdir}/${run}_${cl}_${TAG}" \
    --outSAMtype BAM SortedByCoordinate
done <<< "$SAMPLES"

############### featureCounts ###############
set +u
module unload compiler/GCC 2>/dev/null || true
module add bio/Subread/2.1.1-GCC-13.2.0
set -u

while read -r cl run dir; do
  [ -z "$cl" ] && continue
  sdir="${WORK_DIR}/${dir}"
  bam="${sdir}/${run}_${cl}_${TAG}Aligned.sortedByCoord.out.bam"
  [ -s "$bam" ] || { echo "ERROR: BAM not produced: $bam" >&2; exit 1; }
  echo "=== featureCounts ${run} (${cl}) ==="

  # unique reads only (as in Laura_MANE_personalized.sh)
  featureCounts -T $THREADS -F SAF -p -B -C --countReadPairs -s 0 \
    -a "${REF_DIR}/${cl}_cds_mane.saf" \
    -o "${sdir}/${run}_${cl}_${TAG}_fc.counts.csv" \
    "$bam"

  # multimappers split fractionally (use for HLA allele-pair sums)
  featureCounts -T $THREADS -F SAF -p -B -C --countReadPairs -s 0 -M --fraction \
    -a "${REF_DIR}/${cl}_cds_mane.saf" \
    -o "${sdir}/${run}_${cl}_${TAG}_fcM.counts.csv" \
    "$bam"
done <<< "$SAMPLES"

echo "Done."
exit 0
