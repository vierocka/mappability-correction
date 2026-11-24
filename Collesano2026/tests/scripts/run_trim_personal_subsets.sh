#!/bin/bash
set -euo pipefail
S=${SCRATCH:?set SCRATCH};  I=$S/impact; FQ=$S/hla_test/fq
M=<SELF_NONSELF>/Reference/MANE; PM=$M/personal_MANE_utr; TRIM=<SW>/Trimmomatic
P=<SELF_NONSELF>/Reference/mappability-correction/Collesano2026/scripts
T=10; cd $I
declare -A CELL=([SRR31477874]=HCC1937 [SRR3644141]=HCT116 [SRR364065]=JY)
idx(){ [ -s $1/SA ] || { mkdir -p $1; STAR --runThreadN $T --runMode genomeGenerate --genomeDir $1 --genomeFastaFiles $2 --genomeSAindexNbases 12 --genomeChrBinNbits 11 --outFileNamePrefix $1/ > $1.log; }; }
idx idx_basic $M/MANE.GRCh38.v1.5.ensembl_rna.fna
for c in HCC1937 HCT116 JY; do idx idx_$c $PM/${c}_personal_MANE.fa; done
align(){ # sample arm index saf r1 r2
  local o=$1_$2; [ -s $o.fcM.counts ] && return
  STAR --runMode alignReads --runThreadN $T --genomeDir $3 --readFilesIn $5 $6 --readFilesCommand zcat \
    --alignIntronMax 1 --alignIntronMin 2 --outFilterMultimapNmax 20 --outSAMattributes NH HI AS NM MD \
    --runRNGseed 0 --outSAMtype BAM Unsorted --outFileNamePrefix ${o}_ > /dev/null
  featureCounts -T $T -F SAF -p -B -C --countReadPairs -s 0 -a $4 -o $o.fc.counts ${o}_Aligned.out.bam > /dev/null 2>&1
  featureCounts -T $T -F SAF -p -B -C --countReadPairs -s 0 -M --fraction -a $4 -o $o.fcM.counts ${o}_Aligned.out.bam > /dev/null 2>&1
  rm ${o}_Aligned.out.bam; }
for s in SRR31477874 SRR3644141 SRR364065; do c=${CELL[$s]}; mkdir -p fastqc/$s
  fastqc -q -t 2 --extract -o fastqc/$s $FQ/${s}_1.fastq.gz $FQ/${s}_2.fastq.gz
  # Trimmomatic: PLACEHOLDER = manual standard command (project_Fonly params not available)
  [ -s ${s}_1P.fq.gz ] || java -jar $TRIM/target/trimmomatic-0.40.jar PE -threads $T -phred33 -summary ${s}_trim_summary.txt \
     $FQ/${s}_1.fastq.gz $FQ/${s}_2.fastq.gz ${s}_1P.fq.gz ${s}_1U.fq.gz ${s}_2P.fq.gz ${s}_2U.fq.gz \
     ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> ${s}_trim.log
  mkdir -p fastqc/${s}_trimmed; fastqc -q -t 2 --extract -o fastqc/${s}_trimmed ${s}_1P.fq.gz ${s}_2P.fq.gz
  align $s personal_raw idx_$c $PM/${c}_cds_mane.saf $FQ/${s}_1.fastq.gz $FQ/${s}_2.fastq.gz
  align $s personal_trim idx_$c $PM/${c}_cds_mane.saf ${s}_1P.fq.gz ${s}_2P.fq.gz
  align $s basic_raw idx_basic $M/cds_mane.saf $FQ/${s}_1.fastq.gz $FQ/${s}_2.fastq.gz
  python3 $P/parse_fastqc.py --sample $s --stage raw --fastqc fastqc/$s/*_fastqc/fastqc_data.txt --star-log ${s}_personal_raw_Log.final.out --out ${s}_qc_raw.tsv
  python3 $P/parse_fastqc.py --sample $s --stage trimmed --fastqc fastqc/${s}_trimmed/*_fastqc/fastqc_data.txt --star-log ${s}_personal_trim_Log.final.out --out ${s}_qc_trim.tsv
  echo "done $s"
done
