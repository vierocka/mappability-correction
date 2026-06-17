#!/usr/bin/env bash
# explore_commands.sh — All exploratory bash commands used to design run_hla_mafft.py
# These are one-off queries, not part of the pipeline. Run them individually.

GENCODE=~/Dropbox/Self-Nonself/Reference/GENCODE/gencode.v49.transcripts.fa.gz
MANE=~/Dropbox/Self-Nonself/Reference/MANE/MANE.GRCh38.v1.5.ensembl_rna.fna

# ── 1. Check GENCODE FASTA header format ────────────────────────────────────
zcat "$GENCODE" | grep "^>" | head -3

# ── 2. List all unique HLA gene names in GENCODE v49 ────────────────────────
zcat "$GENCODE" | grep "^>" | awk -F'|' '{print $6}' | grep "^HLA" | sort -u

# ── 3. Longest transcript per HLA pseudogene (any pseudogene biotype) ────────
zcat "$GENCODE" | grep "^>" \
  | awk -F'|' '$6 ~ /^HLA-/ {print $6, $8, $7}' \
  | grep -v "DQB1-AS1\|F-AS1" \
  | awk '{gene=$1; type=$2; len=$3} type ~ /pseudogene/ {print gene, type, len}' \
  | sort -k1,1 -k3,3rn \
  | awk '!seen[$1]++' \
  | sort

# ── 4. Longest protein_coding transcript per functional HLA gene ─────────────
FUNC="HLA-A HLA-B HLA-C HLA-E HLA-F HLA-G HLA-DRB1 HLA-DRB3 HLA-DRB4 HLA-DRB5 HLA-DQA1 HLA-DQB1 HLA-DPA1 HLA-DPB1 HLA-DRA"
zcat "$GENCODE" | grep "^>" \
  | awk -F'|' '$8=="protein_coding" {print $6, $7, $5, $1}' \
  | sort -k1,1 -k2,2rn \
  | awk '!seen[$1]++' \
  | grep -f <(echo "$FUNC" | tr ' ' '\n') \
  | sort

# ── 5. Top HLA-A transcripts by length (to understand isoform size range) ────
zcat "$GENCODE" | grep "^>" \
  | awk -F'|' '$6=="HLA-A" && $8=="protein_coding" {print $7, $5, $1}' \
  | sort -k1,1rn \
  | head -5

# ── 6. HLA MANE Select transcript IDs ───────────────────────────────────────
grep "gene_symbol:HLA-" "$MANE" \
  | grep -oP '>(ENST[0-9.]+).*gene_symbol:(HLA-\S+)' \
  | sed 's/^>//' \
  | awk '{print $NF, $1}' \
  | sed 's/gene_symbol://'

# ── 7. Check mafft version ───────────────────────────────────────────────────
mafft --version

# ── 8. All HLA transcripts with biotype and length (full table) ──────────────
# (output is large; redirect to file)
zcat "$GENCODE" | grep "^>" \
  | awk -F'|' '$6 ~ /^HLA-/ {print $6, $8, $7}' \
  | grep -v "DQB1-AS1\|F-AS1" \
  | sort -k1,1 -k3,3rn \
  > /tmp/hla_transcripts_all.txt

# ── 9. Manual pairwise check: how many transcripts per pseudogene gene ────────
for gene in HLA-H HLA-J HLA-K HLA-L HLA-N HLA-P HLA-S HLA-T HLA-U HLA-V HLA-W HLA-Z \
            HLA-DRB2 HLA-DRB6 HLA-DRB7 HLA-DRB8 HLA-DRB9 HLA-DPA2 HLA-DPB2 HLA-DQB3; do
    n=$(zcat "$GENCODE" | grep "^>" | awk -F'|' -v g="$gene" '$6==g' | wc -l)
    echo "$gene : $n transcripts"
done
