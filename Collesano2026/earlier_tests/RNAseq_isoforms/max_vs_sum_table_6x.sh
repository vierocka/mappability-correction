cd <SELF_NONSELF>/human/RNAseq/Gene_isoform_sample

python3 gene_isoform_sum_max.py \
  --abundance <SELF_NONSELF>/human/RNAseq/HCC1937/kalisto/abundance.tsv \
  --annot <SELF_NONSELF>/human/References/GRCh38.p14/CDS_GRCh38p14_annotNames.csv \
  --out <SELF_NONSELF>/human/RNAseq/HCC1937/kalisto/gene_sum_max_estcounts.tsv \
  --use est_counts
  
python3 gene_isoform_sum_max.py \
  --abundance <SELF_NONSELF>/human/RNAseq/HCC1937/kalisto/abundance.tsv \
  --annot <SELF_NONSELF>/human/References/GRCh38.p14/CDS_GRCh38p14_annotNames.csv \
  --out <SELF_NONSELF>/human/RNAseq/HCC1937/kalisto/gene_sum_max_tpm.tsv \
  --use tpm
  
python3 gene_isoform_sum_max.py \
  --abundance <SELF_NONSELF>/human/RNAseq/JY/abundance.tsv \
  --annot <SELF_NONSELF>/human/References/GRCh38.p14/CDS_GRCh38p14_annotNames.csv \
  --out <SELF_NONSELF>/human/RNAseq/JY/gene_sum_max_estcounts.tsv \
  --use est_counts
  
python3 gene_isoform_sum_max.py \
  --abundance <SELF_NONSELF>/human/RNAseq/JY/abundance.tsv \
  --annot <SELF_NONSELF>/human/References/GRCh38.p14/CDS_GRCh38p14_annotNames.csv \
  --out <SELF_NONSELF>/human/RNAseq/JY/gene_sum_max_tpm.tsv \
  --use tpm
  
python3 gene_isoform_sum_max.py \
  --abundance <SELF_NONSELF>/human/RNAseq/HCT116/abundance.tsv \
  --annot <SELF_NONSELF>/human/References/GRCh38.p14/CDS_GRCh38p14_annotNames.csv \
  --out <SELF_NONSELF>/human/RNAseq/HCT116/gene_sum_max_estcounts.tsv \
  --use est_counts
  
python3 gene_isoform_sum_max.py \
  --abundance <SELF_NONSELF>/human/RNAseq/HCT116/abundance.tsv \
  --annot <SELF_NONSELF>/human/References/GRCh38.p14/CDS_GRCh38p14_annotNames.csv \
  --out <SELF_NONSELF>/human/RNAseq/HCT116/gene_sum_max_tpm.tsv \
  --use tpm
  
  
  
