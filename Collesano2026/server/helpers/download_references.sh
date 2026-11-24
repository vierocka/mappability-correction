#!/bin/bash
# Download and prepare references into refs/ (run once, from Collesano2026/).
#   MANE Select v1.5 (GRCh38): transcripts, GFF, proteins            [Morales 2022]
#   Ensembl 116 GRCh38 primary assembly (194 seqs; funique only)       [Ensembl]
#   IPD-IMGT/HLA 3.63.1 hla_nuc.fasta                                  [Robinson 2019]
#   Trimmomatic TruSeq3-PE adapters                                     [Bolger 2014]
#   arcasHLA reference rebuilt from the same IMGT release               [Orenbuch 2020]
set -euo pipefail
R=${1:-refs}; mkdir -p $R/{MANE,Ensembl,IMGTHLA,adapters}
MANE=https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5
for f in MANE.GRCh38.v1.5.ensembl_rna.fna MANE.GRCh38.v1.5.ensembl_genomic.gff MANE.GRCh38.v1.5.ensembl_protein.faa; do
  [ -s $R/MANE/$f ] || { wget -q -O - $MANE/$f.gz | gunzip -c > $R/MANE/$f; }
done
G=Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
[ -s $R/Ensembl/$G ] || wget -q -O - https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/$G.gz | gunzip -c > $R/Ensembl/$G
n=$(grep -c '>' $R/Ensembl/$G); [ "$n" -eq 194 ] || { echo "expected 194 sequences, got $n"; exit 1; }
IMGT=https://raw.githubusercontent.com/ANHIG/IMGTHLA/3631
[ -s $R/IMGTHLA/hla_nuc.fasta ] || wget -q -O $R/IMGTHLA/hla_nuc.fasta $IMGT/hla_nuc.fasta
[ -s $R/IMGTHLA/hla.dat ] || { wget -q -O $R/IMGTHLA/hla.dat.zip $IMGT/hla.dat.zip && (cd $R/IMGTHLA && unzip -o -q hla.dat.zip); }
[ -s $R/adapters/TruSeq3-PE.fa ] || wget -q -O $R/adapters/TruSeq3-PE.fa \
  https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa
python3 "$(dirname "$0")/../../scripts/make_mane_saf.py" $R/MANE/MANE.GRCh38.v1.5.ensembl_rna.fna $R/MANE/cds_mane.saf
echo "References ready in $R."
echo "arcasHLA reference (once, in the arcashla env):"
echo "  build a local git snapshot of IMGTHLA 3.63.1 in \$(dirname \$(which arcasHLA))/../share/arcas-hla-0.6.0-*/dat/IMGTHLA"
echo "  (hla.dat + wmda/hla_nom_g.txt + wmda/hla_nom_p.txt), then: arcasHLA reference --rebuild"
