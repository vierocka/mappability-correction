# source ~/miniconda3/bin/activate arcasHLA
# arcasHLA reference --rebuild

BASEDIR="<SELF_NONSELF>/arcasHLA"

cd $BASEDIR

# Sample directories and BAM prefixes
SAMPLES=(
    "SRR364065_JY"
    "SRR31477875_HCC1937"
    "SRR31477874_HCC1937"
    "SRR925709_HCC1937"
    "SRR3644141_HCC116"
    "SRR3644142_HCC116"
    "SRR3644143_HCC116"
    "SRR3644144_HCC116"
)

for SAMPLE in "${SAMPLES[@]}"; do

    # Parse directory and prefix
    PREFIX="${SAMPLE#*/}"
    WORKDIR="${BASEDIR}"

    echo "=============================================="
    echo "Processing: ${PREFIX}"
    echo "Directory:  ${WORKDIR}"
    echo "=============================================="

# 2.) GENOTYPE - identify HLA alleles
    echo "[$(date)] Genotyping..."
    arcasHLA genotype \
        ${PREFIX}_gencodeAligned.sortedByCoord.out.extracted.1.fq.gz \
        ${PREFIX}_gencodeAligned.sortedByCoord.out.extracted.2.fq.gz \
        -g A,B,C \
	-t 8 \
	--max_iterations 500 \
	--min_count 5 \
        -o ./

    # 3.) QUANTIFY - estimate allele abundance
    echo "[$(date)] Quantifying allele expression..."
    arcasHLA quant \
        ${PREFIX}_gencodeAligned.sortedByCoord.out.extracted.1.fq.gz \
        ${PREFIX}_gencodeAligned.sortedByCoord.out.extracted.2.fq.gz \
       # --genotype ${PREFIX}_gencodeAligned.sortedByCoord.out.genotype.json \
        -t 8 \
        -o ./

    echo "[$(date)] Done with ${PREFIX}"
    echo ""

done

echo "=============================================="
echo "All samples processed!"
echo "=============================================="

# Summary: collect all genotypes
echo ""
echo "Genotype summary:"
for SAMPLE in "${SAMPLES[@]}"; do
    PREFIX="${SAMPLE#*/}"
    JSON="${BASEDIR}/${PREFIX}_gencodeAligned.sortedByCoord.out.genotype.json"
    if [ -f "${JSON}" ]; then
        echo "${PREFIX}: $(cat ${JSON})"
    else
        echo "${PREFIX}: GENOTYPE NOT FOUND"
    fi
done
