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
| `--outFilterMultimapNmax` | 20 (was 40 before 2026-07-08) | Matches nf-core/rnaseq's ENCODE3-standard `star_rsem` preset (see `../../Ref/recommendations_nf.txt`). Changed across all scripts, old and `_nf` alike — see `roadmap.csv`'s `outFilterMultimapNmax` column. Retains multimappers in BAM for diagnostic purposes up to that cap. |
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

## nf reference variant (`_nf` scripts)

42 scripts re-run three arms of the pipeline against the Ensembl GRCh38.116
primary assembly + matching GTF recommended by nf-core/rnaseq, instead of
the NCBI RefSeq GRCh38.p14 reference used above — to check how much of the
uniqueness-factor signal is reference-choice-dependent:

- `v1_a-d_nf`, `v3_a-d_nf`, `v4_a-d_nf`, `v6_a-b_nf`, `v8_a_nf` (15) — genomic reference, MANE_Select scope. See below for what differs from the non-`_nf` scripts.
- `v7_a-j_nf`, `v8_b_nf`, `v8_c_nf` (12) — genome-wide dedup-CDS reference, added once it became clear a fair comparison against the old `v7`/`v8_b` results needs the same genome-wide (not MANE-restricted) CDS scope — see the dedicated subsection below.
- `v2_a-o_nf` (15) — MANE-restricted dedup-CDS reference, a **third, separate** CDS reference from the genome-wide one above — see its own dedicated subsection below for why.

`v5` (MANE transcriptome as a standalone self-map reference target) still
doesn't have an `_nf` counterpart — `MANE_RNA_nf.fna` already exists (built
by Step 0c, shared with `v8_c_nf`), so this would be a comparatively small
addition if wanted later; see `run_all_mappability_nf.sh`'s header for the
current, authoritative list.

**What differs from the scripts above (genomic-reference arm, `v1/v3/v4/v6/v8_a` family):**
- Reference: `Ref/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa` + `Ref/Homo_sapiens.GRCh38.116.gtf` — lives inside `../../Ref/` itself, alongside the original NCBI RefSeq files, distinguished by filename (`Homo_sapiens.GRCh38.*` vs. `GCF_*`/`MANE.*`) rather than a separate directory. Nothing here ever touches the original NCBI files, indices, or non-`_nf` scripts.
- Exon/transcript coordinates come from the Ensembl GTF filtered to `tag "MANE_Select"` transcripts, not the separate MANE GFF — the GTF and genome FASTA already share chromosome naming (`1`, `2`, …, `X`, `Y`, `MT`), so no chr-name translation is needed internally (unlike v1_a.sh's Step 0). `chr_name_map_nf.tsv` is still produced, but only as a cross-reference artifact for comparing old vs. new results by chromosome (see `../results/README.md`).
- MANE_RNA-sourced scripts need spliced transcript sequences for the new reference; since nf-core recommends deriving these rather than downloading a separate transcriptome FASTA (see `../../Ref/recommendations_nf.txt`), Step 0c splices MANE_Select transcripts from the GTF + genome into `Ref/MANE_RNA_nf.fna`, built once and shared by every MANE-RNA-sourced `_nf` script.
- Each `_nf` script is standalone (own `module purge`/`module load`, own `ulimit -n 65536`) since there is no `_nf` equivalent of the old convention where individual scripts relied on a driver to load modules once.
- Generated by `generate_scripts_nf.py` (all except `v1_a_nf.sh`, hand-written as the reference pattern the generator replicates).

### Genome-wide dedup-CDS arm (`v7_a-j_nf`, `v8_b_nf`, `v8_c_nf`)

Unlike every other `_nf` script, this arm deliberately is **not** restricted
to MANE_Select — the old `v2`/`v7`/`v8_b` dedup-CDS reference covers every
protein-coding transcript genome-wide (93,088 sequences after dedup), and
narrowing that to MANE_Select for the nf version would make it a different,
smaller comparison rather than a fair like-for-like against the old results.

- `Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna` is built by `v7_a_nf`'s Step 0d: parses `Ref/Homo_sapiens.GRCh38.116.gtf` for `CDS` features on `transcript_biotype "protein_coding"` transcripts (not MANE-filtered), splices multi-exon CDS per transcript via `pysam.FastaFile` (not per-exon `samtools` subprocess calls — this is a much larger extraction than the MANE-only one), reverse-complements minus-strand transcripts, and writes `Ref/Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_raw.fna` with plain versioned-transcript-ID headers (e.g. `ENST00000256078.10`) — no `|` characters, unlike NCBI's `lcl|...` CDS IDs. That choice sidesteps the bug `mappability_v7_fix_uf.sh` had to patch post-hoc in the old pipeline (grouping all reads under `cds_id='lcl'` because the ID itself contained `|`).
- Deduplication (`seqkit rmdup -s`) has no confirmed module on this cluster — Step 0d checks for `seqkit` on `PATH` after extraction; if it's missing, it stops with instructions to `scp` `Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_raw.fna` to a machine with seqkit installed (any recent version — `rmdup -s` has been stable for years), run the dedup there, and `scp Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna` back to the same path in `Ref/`. Re-running then just skips straight past Step 0d.
- `v7_a_nf` also builds `STAR_index_dedup_cds_nf`; `v8_c_nf` builds `STAR_index_MANE_nf` (from the already-existing `MANE_RNA_nf.fna`).
- `v7_a-f_nf`/`v7_i_nf`/`v7_j_nf` align CDS-derived reads to `STAR_index_genome_nf`; `v7_g_nf`/`v7_h_nf` align to `STAR_index_dedup_cds_nf` itself (self-map — residual within-CDS redundancy after dedup).
- `v8_b_nf`/`v8_c_nf` reuse `v8_a_nf`'s simulated reads (checked explicitly; they abort with a clear error if `v8_a_nf` hasn't been run yet) and align to `STAR_index_dedup_cds_nf` / `STAR_index_MANE_nf` respectively, with `--alignIntronMax 1` (matching the old `v8_b`/`v8_c`).
- Output columns match the old `cds_uniqueness_factors_*.tsv` schema (`cds_id`, not `transcript_id`) for `v7_*_nf`; `v8_b_nf`/`v8_c_nf` use the `transcript_id` schema, matching old `v8_b`/`v8_c`.

### MANE-restricted dedup-CDS arm (`v2_a-o_nf`)

The old `v2` scripts align MANE-sourced reads (exon-flank / MANE RNA) against
the **same genome-wide** `dedup_cds.fna` used by `v7` — confirmed directly:
`v2_a.sh`/`v2_c.sh` both point at `STAR_index_dedup_cds`. So a literal `v2_nf`
mirror would just reuse `Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna`/`STAR_index_dedup_cds_nf` from the
arm above. Instead, `v2_nf` deliberately uses a **third, separate** CDS
reference restricted to MANE_Select transcripts only — a different, smaller
comparison than either old `v2` or genome-wide `v7_nf`/`v8_b_nf`.

- `Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna` is built by `v2_a_nf`'s Step 0d: same splicing approach as `Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna` (via `pysam.FastaFile`, plain versioned-transcript-ID headers), but parses `CDS` features only on `tag "MANE_Select"` transcripts (~19,400, reusing the same 2-pass tag-scan as `MANE_RNA_nf.fna`'s Step 0c) instead of every `transcript_biotype "protein_coding"` transcript genome-wide.
- Same resumable seqkit fallback as the genome-wide build: extraction always runs (pure Python + pysam), dedup only runs if `seqkit` is on `PATH`, otherwise it stops with `scp`/dedup/`scp`-back instructions for `Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane_raw.fna` → `Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna`.
- `v2_a_nf` also builds `STAR_index_dedup_cds_mane_nf`.
- Read simulation is identical to the genomic-reference arm's exon-flank/MANE-RNA logic (Step 1 is reused unchanged) — only the STAR alignment target differs (`STAR_index_dedup_cds_mane_nf` instead of `STAR_index_genome_nf`).
- `v2_d_nf` sets `--alignIntronMax 1`, matching old `v2_d`.
- Naming mirrors the old `v2` scripts' idiosyncrasies exactly (plus `_nf`) — e.g. `v2_c_nf`'s bare `dedup_cds_nf` output dir and `v2_e_nf`'s missing `L100bp` tag, both inherited from the original hand-authored `v2_c.sh`/`v2_e.sh` naming.

**Usage:** run all 42 via the driver, `bash run_all_mappability_nf.sh` (submit with `sbatch` on the cluster) — order: `v1_a_nf`/`v1_c_nf` (shared `STAR_index_genome_nf` + `MANE_RNA_nf.fna`), rest of the genomic arm, `v7_a_nf` (shared `Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna` + `STAR_index_dedup_cds_nf`), rest of `v7`, `v8_b_nf`/`v8_c_nf`, then `v2_a_nf` (shared `Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna` + `STAR_index_dedup_cds_mane_nf`) and the rest of `v2`. Skips any script whose output TSV already exists so a partial run can be resubmitted safely — including safely resuming after a manual seqkit round-trip for either CDS reference's Step 0d.

Annotation cross-checks (`check_id_conversion_nf.py`, `check_exon_positions_nf.py`, one directory up) are independent of these scripts — they only need the GTF, not the genome FASTA or STAR — and are worth running first to sanity-check the MANE_Select transcript set before committing cluster time to the full sweep. `extract_dedup.log` (one directory up) records CDS extraction/dedup counts for both the genome-wide and MANE-restricted builds, old vs. new reference.
