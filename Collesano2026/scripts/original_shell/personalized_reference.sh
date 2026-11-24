# 1. Remove all HLA class I from deduplicated FASTA
cd <SELF_NONSELF>/human/References/GRCh38.p14
seqkit grep -v -r -p "gene=HLA-[ABC]" GCF_000001405.40_GRCh38.p14_cds_from_genomic.dedup.fna > no_HLA_classI.fna
    
# 2. Get CDS for confirmed alleles from IMGT/HLA
# Extract from the nucleotide FASTA
cd <LOCAL>/IMGTHLA-Latest

# JY (homozygous at all loci)
awk 'BEGIN{RS=">"} /A\*02:01:01:01 / {print ">"$0}' fasta/A_nuc.fasta > <SELF_NONSELF>/Reference/GRC38p14/JY_alleles.fa
awk 'BEGIN{RS=">"} /B\*07:02:01:01 / {print ">"$0}' fasta/B_nuc.fasta >> <SELF_NONSELF>/Reference/GRC38p14/JY_alleles.fa
awk 'BEGIN{RS=">"} /C\*07:02:01:03 / {print ">"$0}' fasta/C_nuc.fasta >> <SELF_NONSELF>/Reference/GRC38p14/JY_alleles.fa

# HCT116 (heterozygous)
awk 'BEGIN{RS=">"} /A\*01:01:01:01 / {print ">"$0}' fasta/A_nuc.fasta > <SELF_NONSELF>/Reference/GRC38p14/HCT116_alleles.fa
awk 'BEGIN{RS=">"} /A\*02:01:01:01 / {print ">"$0}' fasta/A_nuc.fasta >> <SELF_NONSELF>/Reference/GRC38p14/HCT116_alleles.fa
awk 'BEGIN{RS=">"} /B\*18:01:01:01 / {print ">"$0}' fasta/B_nuc.fasta >> <SELF_NONSELF>/Reference/GRC38p14/HCT116_alleles.fa
awk 'BEGIN{RS=">"} /B\*45:01:01:01 / {print ">"$0}' fasta/B_nuc.fasta >> <SELF_NONSELF>/Reference/GRC38p14/HCT116_alleles.fa
awk 'BEGIN{RS=">"} /C\*05:01:01:01 / {print ">"$0}' fasta/C_nuc.fasta >> <SELF_NONSELF>/Reference/GRC38p14/HCT116_alleles.fa
awk 'BEGIN{RS=">"} /C\*07:01:01:01 / {print ">"$0}' fasta/C_nuc.fasta >> <SELF_NONSELF>/Reference/GRC38p14/HCT116_alleles.fa

# HCC1937 (heterozygous)
awk 'BEGIN{RS=">"} /A\*23:01:01:01 / {print ">"$0}' fasta/A_nuc.fasta > <SELF_NONSELF>/Reference/GRC38p14/HCC1937_alleles.fa
awk 'BEGIN{RS=">"} /A\*24:02:01:01 / {print ">"$0}' fasta/A_nuc.fasta >> <SELF_NONSELF>/Reference/GRC38p14/HCC1937_alleles.fa
awk 'BEGIN{RS=">"} /B\*07:02:01:01 / {print ">"$0}' fasta/B_nuc.fasta >> <SELF_NONSELF>/Reference/GRC38p14/HCC1937_alleles.fa
awk 'BEGIN{RS=">"} /B\*40:01:01/ {print ">"$0}' fasta/B_nuc.fasta >> <SELF_NONSELF>/Reference/GRC38p14/HCC1937_alleles.fa
awk 'BEGIN{RS=">"} /C\*03:04:01:01 / {print ">"$0}' fasta/C_nuc.fasta >> <SELF_NONSELF>/Reference/GRC38p14/HCC1937_alleles.fa
awk 'BEGIN{RS=">"} /C\*07:02:01:01 / {print ">"$0}' fasta/C_nuc.fasta >> <SELF_NONSELF>/Reference/GRC38p14/HCC1937_alleles.fa

### Verify
### grep -c ">" <SELF_NONSELF>/Reference/GRC38p14/JY_alleles.fa      # should be 3
### grep -c ">" <SELF_NONSELF>/Reference/GRC38p14/HCT116_alleles.fa   # should be 6

# 3. Append correct alleles per cell line
cat <SELF_NONSELF>/human/References/GRCh38.p14/no_HLA_classI.fna <SELF_NONSELF>/Reference/GRC38p14/JY_alleles.fa > personalized_JY.fna
cat <SELF_NONSELF>/human/References/GRCh38.p14/no_HLA_classI.fna <SELF_NONSELF>/Reference/GRC38p14/HCT116_alleles.fa > personalized_HCT116.fna
cat <SELF_NONSELF>/human/References/GRCh38.p14/no_HLA_classI.fna <SELF_NONSELF>/Reference/GRC38p14/HCC1937_alleles.fa > personalized_HCC1937.fna

# 4. Index + full remap
### hpc - copy the personalized references
# kallisto index -i personalized_JY.idx personalized_JY.fna
# kallisto quant -i personalized_JY.idx -t 24 -b 100 -o JY_personalized SRR364065_1.fastq SRR364065_2.fastq




    


