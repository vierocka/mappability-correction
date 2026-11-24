#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --mem=192gb 
#SBATCH --time=4:00:00 
#SBATCH --output=HCC1937_kallisto.out 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>
#SBATCH --account=<group>

HumanRef="<HPC_WORKDIR>/GRCh38.p14"
ColonCancerSample="<HPC_WORKDIR>/Laura/"

cd <HPC_WORKDIR>/GRCh38.p14

############ CHANGE OF COMPILER ########
module unload compiler/GCC
module add bio/kallisto/0.50.1-gompi-2022b

### colon cancer: HCC1937 - SRR925709 - done
### breast cancer HCT116 - SRR15050724
cd <HPC_WORKDIR>/Laura/SRR15050724
mkdir kallisto_quant
kallisto quant -i $HumanRef/kallisto_index_GRCh38p14.idx -o kallisto_quant -b 100 -t 24 SRR15050724_1P.fq SRR15050724_2P.fq

### JY - SRR364065 - JY has many unconventional isoforms
cd <HPC_WORKDIR>/Laura/SRR364065
mkdir kallisto_quant
kallisto quant -i $HumanRef/kallisto_index_GRCh38p14.idx -o kallisto_quant -b 100 -t 24 SRR364065_1P.fq SRR364065_2P.fq

### breast cancer
# SRR925689 - 21MT1 - metastatic - strong “viral-mimic/ISG” immune-mimic signatures
# SRR925690 - 21MT2 - metastatic - strong “viral-mimic/ISG” immune-mimic signatures
# SRR925691 - 21NT1 - Primary non-tumorigenic mammary epithelial cell line, patient 21 (preneoplastic / premalignant epithelium); “normal adjacent” epithelium
# SRR925692 - 21 PT- primary tumor

cd <HPC_WORKDIR>/Laura/SRR925689
mkdir kallisto_quant
kallisto quant -i $HumanRef/kallisto_index_GRCh38p14.idx -o kallisto_quant -b 100 -t 24 SRR925689_1P.fq SRR925689_2P.fq

cd <HPC_WORKDIR>/Laura/SRR925690
mkdir kallisto_quant
kallisto quant -i $HumanRef/kallisto_index_GRCh38p14.idx -o kallisto_quant -b 100 -t 24 SRR925690_1P.fq SRR925690_2P.fq

cd <HPC_WORKDIR>/Laura/SRR925691
mkdir kallisto_quant
kallisto quant -i $HumanRef/kallisto_index_GRCh38p14.idx -o kallisto_quant -b 100 -t 24 SRR925691_1P.fq SRR925691_2P.fq

cd <HPC_WORKDIR>/Laura/SRR925692
mkdir kallisto_quant
kallisto quant -i $HumanRef/kallisto_index_GRCh38p14.idx -o kallisto_quant -b 100 -t 24 SRR925692_1P.fq SRR925692_2P.fq

exit 0

