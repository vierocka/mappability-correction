# Map_scripts — simulation and alignment scripts

Scripts for simulating reads from reference sequences and computing per-transcript uniqueness factors with STAR. See `../roadmap.csv` for the complete parameter matrix and expected output paths for all 56 variants.

## Versioning scheme

| Version family | Read source | Reference targets | Key variable |
|---|---|---|---|
| v1 | Genomic exon ± 50 bp flanks or MANE RNA | Full genome | Mode (PE/SE), read source |
| v2 | Genomic exon ± 50 bp flanks or MANE RNA | Deduplicated CDS (protein-coding only) | Mode, read source, alignIntronMax |
| v3 | Genomic exon ± 50 bp flanks | Full genome | Read length (75 bp) |
| v4 | Genomic exon ± 50 bp flanks or MANE RNA | Full genome | Read length (150 bp) |
| v5 | MANE RNA or genomic exon flanks | MANE transcriptome | Mode, read source, read length |
| v6 | MANE RNA | Full genome | Read length (200 bp) |
| v7 | Deduplicated CDS sequences (protein-coding only; no reads simulated for non-coding genes) | Full genome or dedup CDS self-map | Read length, mode |
| v8 | MANE RNA, proper PE 300 bp insert | Full genome / dedup CDS / MANE | Reference (a=genome, b=CDS, c=MANE) |

## Selected parameter comparisons

Key pairwise comparisons and what they isolate:

| Comparison | What changes | What is isolated |
|---|---|---|
| v1_a vs v1_b | PE → SE | Concordance filter effect |
| v1_b vs v1_c | Genomic flanks → spliced RNA reads | Intronic flank contribution |
| v1_b vs v2_b | Full genome → dedup CDS | Reference competing-sequence content |
| v2_c vs v2_d | No alignIntronMax vs alignIntronMax 1 | STAR splice-search artefact on CDS reference |
| v1_a / v3_a / v4_a | 100 → 75 → 150 bp | Read length effect on UF |
| v5 vs v1/v2 | MANE transcriptome ref | Self-mapping completeness |
| v7_a vs v7_g | Full genome vs dedup CDS self-map | Residual within-CDS redundancy |
| v4_d vs v8_a | Fake PE (R2=revcomp R1) vs proper PE | Insert geometry effect |

## v8 — proper paired-end geometry

v1–v7 simulate pairs where R2 = revcomp(R1) at the same window (insert size = 0), which is an ideal model read. v8 simulates a realistic 300 bp fragment: R1 is the first 150 bp, R2 is the reverse complement of the last 150 bp. This matches real Illumina PE 150 bp library structure and activates STAR's proper-pair concordance checking.

v8_b and v8_c **reuse the simulated reads from v8_a** — run v8_a first. The shared read directory is `../results/v8_simreads_MANE_RNA_PE_L150bp_proper/`.

STAR parameters differ by reference type:

| Script | Reference | alignIntronMax | Reason |
|---|---|---|---|
| v8_a | Full genome | (not set) | Full genome needs splice-spanning alignment |
| v8_b | Dedup CDS | 1 | CDS sequences are already unspliced exonic sequence |
| v8_c | MANE RNA | 1 | Spliced transcripts — no introns to span |

Index build parameters for small references (dedup CDS, MANE RNA):
- `--genomeSAindexNbases 12`
- `--genomeChrBinNbits 11`

## STAR alignment parameters (all scripts)

| Parameter | Value | Reason |
|---|---|---|
| `--alignEndsType` | Local | Reads may not cover the full fragment end; local alignment avoids soft-clip penalties inflating MAPQ |
| `--outFilterMultimapNmax` | 100 | Retain multimappers in BAM for diagnostic purposes |
| `--alignMatesGapMax` | 0 | No gap limit between mates (0 = unlimited); transcript lengths vary |
| `--outSAMattributes` | NH HI AS NM | NH (total alignments) needed to identify unique mappers at MAPQ=255 |

Uniqueness factor = MAPQ-255 R1 reads / total simulated R1 reads per transcript.

## Dependencies

- [STAR](https://github.com/alexdobin/STAR) ≥ 2.7
- [samtools](https://www.htslib.org/) ≥ 1.15
- Python ≥ 3.9 with `pysam` and `pandas`
- [seqkit](https://bioinf.shenwei.me/seqkit/) — for CDS deduplication only (run once, not per script)

## Reference files

The `../../Ref/` directory is not tracked by git. Populate it as follows:

```bash
# Full genome (GRCh38.p14) — download from NCBI FTP
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/
# File: GCF_000001405.40_GRCh38.p14_genomic.fna (~3.2 GB)

# CDS sequences — download, then deduplicate
seqkit rmdup -s -i GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna \
    -o Ref/GCF_000001405.40_GRCh38.p14_cds_from_genomic.dedup.fna \
    --threads 12

# MANE Select v1.5 RNA sequences and GFF — download from NCBI MANE FTP
# https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/
# Files: MANE.GRCh38.v1.5.ensembl_rna.fna, MANE.GRCh38.v1.5.ensembl_genomic.gff
```

A pre-built STAR index for the full genome (with MANE GFF annotation) is required before running v1–v4, v6, v7, and v8_a. The dedup CDS and MANE RNA indices are built by the scripts on first run.

See `../../ref_files_needed.csv` for the complete file inventory.

## Usage

Edit the path variables at the top of each script, then run sequentially:

```bash
bash mappability_correction_v1_a.sh
```

For v8, run in order:
```bash
bash mappability_correction_v8_a.sh   # generates reads + aligns to full genome
bash mappability_correction_v8_b.sh   # reuses reads, aligns to dedup CDS
bash mappability_correction_v8_c.sh   # reuses reads, aligns to MANE RNA
```

Each script writes results to `../results/<subdir>/` and produces a `transcript_uniqueness_factors_*.tsv` with one row per transcript.

For HPC submission, adapt `kallisto_mappability_crossval.sh` as a template — it shows the module load commands and SLURM resource requests used on the original cluster.
