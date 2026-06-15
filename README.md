# mappability-correction
RNA-seq analysis assumes that annotated transcripts can be quantified accurately from short reads. This project tests that assumption. Using systematic simulations across references, read lengths and quantification strategies, I estimate transcript-specific uniqueness factors and identify genes that are reliably quantifiable, partially recoverable, or fundamentally unresolved from short-read RNA-seq.

Long-read sequencing can resolve many transcripts that remain ambiguous in short-read RNA-seq. However, its current cost limits routine application in large transcriptomic studies. By identifying problematic loci, this project may help guide targeted long-read sequencing while retaining cost-effective Illumina sequencing for the remainder of the transcriptome.

Many conclusions in immuno-oncology and transcriptomics are based on RNA-seq measurements generated from short-read data. By identifying loci that are systematically misquantified due to pseudogenes, paralogs or annotation complexity, this project may improve both the interpretation of historical datasets and the design of future hybrid short-read/long-read studies.

Despite the rapid adoption of long-read technologies, most large-scale transcriptomic resources—including TCGA, GTEx, TARGET, PCAWG, CCLE, DepMap, ENCODE and the Human Cell Atlas—were generated using short-read RNA-seq and will remain a cornerstone of biomedical research for years to come.

## What this project investigates

Standard RNA-seq pipelines produce opposing systematic biases at loci with complex annotation geometry — dense pseudogene neighbourhoods, gene family duplications, and repetitive elements:

- **Signal loss** — alignment-based counting (STAR + featureCounts) discards multimapping reads, underestimating expression where competing sequences are present in the reference.
- **Signal inflation** — pseudoaligners (Kallisto, Salmon) collapse reads from unannotated or pseudogene sources onto the nearest reference transcript, overestimating expression.
- **Reference dependency** — the direction and magnitude of bias depend on which reference is used (full genome, CDS-only, or spliced transcriptome), not just on the tool.
- **RNA–protein discrepancy** — both biases contribute to the poor mRNA–protein correlation observed at affected loci in proteogenomics datasets.

The project computes per-transcript **uniqueness factors** across a systematic parameter sweep (reference × read length × mode × read source) and uses them to classify every MANE Select transcript into one of four reliability classes: always reliable, apply correction, tool-specific, or unquantifiable from short reads.

For rationale, affected gene categories, and the classification framework, see `mappability/README.md`.

## What is tested

Reads simulated from all MANE Select transcripts (script families v1–v6, v8) and from deduplicated CDS sequences (v7) are aligned back to three references with STAR. The fraction of reads recovering MAPQ = 255 is the per-transcript uniqueness factor. Kallisto cross-validates the inflation direction; featureCounts cross-validates count-level distortion.

The full parameter sweep (56 variants) is listed in `mappability/roadmap.csv`.

## Repository structure

| Folder / file | Contents |
|---|---|
| `mappability/` | Simulation scripts, results, and documentation; roadmap in `roadmap.csv` |
| `HLA_MAFFT/` | MAFFT multiple alignments of HLA/APM sequences across three reference sets; sequence identity figures |
| `Figure_APMlocus/` | Annotation map of the APM/HLA region on GRCh38.p14 and the script to generate it |
| `Ref/` | Reference files — not tracked; see `ref_files_needed.csv` for download sources |
| `comparison_scheme.md` | Comparison design flowchart: what is compared, why, and the escalation path from RNA-seq to immunopeptidomics |
| `ref_files_needed.csv` | Reference file inventory with download sources |
| `deduplication_note.txt` | CDS deduplication command and statistics |
| `versions/` | Tool version log |

## Status

| Module | Status |
|---|---|
| APM locus annotation figure | Done |
| HLA sequence identity (MAFFT) | Done |
| STAR simulation scripts v1–v8 | In progress |
| featureCounts cross-validation | Planned |
| Kallisto cross-validation | Pending |
| Benchmarking and correction factor tables | Pending |
