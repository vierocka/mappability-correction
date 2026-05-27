# mappability-correction

In silico mappability correction for RNA-seq — MANE Select transcripts.

Simulates reads from all 19,373 MANE Select v1.5 transcripts, aligns them back to a reference with STAR, and computes a per-transcript **uniqueness factor**: the fraction of simulated positions that map back uniquely (MAPQ = 255). These factors are used to correct TPM values for transcripts in low-mappability regions — segmental duplications, large gene families (olfactory receptors, HLA), and pseudogene-dense loci.

At 100 bp read length against the full GRCh38.p14 genome, 5.6% of MANE Select transcripts (1,079 / 19,373) have a uniqueness factor of zero — every simulated read multimaps.

## Design

The scripts systematically vary four parameters to assess their effect on uniqueness factor estimates:

| Script | Read len | Reference | Read source | Flanks | Mode | alignIntronMax 1 |
|--------|----------|-----------|-------------|--------|------|-----------------|
| v1_a | 100 bp | Full genome + alt loci | Genomic exon ± 50 bp | Yes | PE | No |
| v1_b | 100 bp | Full genome + alt loci | Genomic exon ± 50 bp | Yes | SE | No |
| v1_c | 100 bp | Full genome + alt loci | MANE spliced RNA | No | SE | No |
| v2_a | 100 bp | Dedup CDS | Genomic exon ± 50 bp | Yes | PE | No |
| v2_b | 100 bp | Dedup CDS | Genomic exon ± 50 bp | Yes | SE | No |
| v2_c | 100 bp | Dedup CDS | MANE spliced RNA | No | SE | No |
| v2_d | 100 bp | Dedup CDS | MANE spliced RNA | No | SE | Yes |
| v3_a | 75 bp | Full genome + alt loci | Genomic exon ± 50 bp | Yes | PE | No |
| v4_a | 150 bp | Full genome + alt loci | Genomic exon ± 50 bp | Yes | PE | No |

**Key pairwise comparisons:**
- `v1_a → v1_b` — PE vs SE
- `v1_b → v1_c` — intronic flanks vs spliced RNA source
- `v1_b → v2_b` — full genome vs deduplicated CDS (cleanest reference comparison)
- `v1_c → v2_c` — same, using spliced RNA reads
- `v2_c → v2_d` — effect of disabling STAR de novo splice search on CDS reference
- `v1_a / v3_a / v4_a` — read length effect (75 / 100 / 150 bp)

**References:**
- *Full genome*: GRCh38.p14 (`GCF_000001405.40_GRCh38.p14_genomic.fna`) — includes primary chromosomes and alt loci
- *Dedup CDS*: CDS sequences from GRCh38.p14 (`GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna`), deduplicated with `seqkit rmdup -s` (53,249 exact duplicates removed, 93,088 unique sequences retained)
- *MANE RNA*: `MANE.GRCh38.v1.5.ensembl_rna.fna` — spliced transcript sequences for all MANE Select transcripts

## Dependencies

- [STAR](https://github.com/alexdobin/STAR) ≥ 2.7
- [samtools](https://www.htslib.org/) ≥ 1.15
- [seqkit](https://bioinf.shenwei.me/seqkit/) (for CDS deduplication only)
- Python ≥ 3.9 with `pysam` and `pandas`

## Reference files

The `Ref/` directory is not tracked by git (files exceed GitHub's size limits). Populate it as follows:

```bash
# Full genome — download from NCBI
# GCF_000001405.40_GRCh38.p14_genomic.fna  (~3.2 GB)
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/

# CDS sequences — download from NCBI, then deduplicate
# GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna
seqkit rmdup -s -i GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna \
    -o Ref/GCF_000001405.40_GRCh38.p14_cds_from_genomic.dedup.fna \
    --threads 12

# MANE Select v1.5 RNA sequences and genomic GFF — download from NCBI MANE FTP
# MANE.GRCh38.v1.5.ensembl_rna.fna
# MANE.GRCh38.v1.5.ensembl_genomic.gff
# https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/
```

You also need a pre-built STAR genome index for the full GRCh38.p14 genome (splice-aware, with the MANE GFF annotation). The CDS STAR index is built automatically by the v2 scripts on first run.

## Usage

Edit the path variables at the top of each script to match your environment, then run:

```bash
bash mappability_correction_v1_a.sh
```

Each script writes results to a dedicated output directory and produces a `transcript_uniqueness_factors_*.tsv` with one row per MANE transcript.
