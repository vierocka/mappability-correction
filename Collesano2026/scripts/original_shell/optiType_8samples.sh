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

cd <SELF_NONSELF>/arcasHLA

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing: ${SAMPLE}"
    OptiTypePipeline.py \
        -i ${SAMPLE}_gencodeAligned.sortedByCoord.out.extracted.1.fq.gz \
           ${SAMPLE}_gencodeAligned.sortedByCoord.out.extracted.2.fq.gz \
        --rna \
        -v \
        -o ./optitype_${SAMPLE}
    echo "Done: ${SAMPLE}"
done
