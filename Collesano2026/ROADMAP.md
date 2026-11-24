# Architecture of `Collesano2026/`

What is in this folder and why. Scope: RNA-seq → corrected TPM → 9-mer source repertoire.
The kinetic model and the mass-spectrometry integration live in a companion repository.

## Folders

| folder | what it is | why it is here |
|---|---|---|
| `server/` | Snakemake workflow: FASTQ → QC → HLA typing → personal MANE → STAR + featureCounts → corrected TPM | the cluster half of the pipeline; see [HOW2RUN.md](HOW2RUN.md) |
| `local/` | Snakemake workflow: corrected TPM + personal proteome → 9-mer net sums | the laptop half; needs `server` output first |
| `scripts/` | the Python steps both workflows call | correction, 9-mer tables, reports, tests |
| `config/` | `server.yaml`, `local.yaml`, `samples.tsv` | all paths and thresholds; one row per sequencing run |
| `envs/` | conda environments per rule group | pinned tool versions (`--use-conda`) |
| `resources/` | precomputed funique tables, non-coding transcript list | so the 40 GB STAR index step can be skipped |
| `tests/` | impact tests: trimming, personal MANE, filters, inflation tail, HLA typing, library composition | each answers "does this choice change the result" |
| `pseudogene_problem/` | the four causes of low uniqueness and which one inflates | why funique exists at all — [the core argument](pseudogene_problem/README.md) |
| `excluded_transcripts/` | the 117 excluded MANE transcripts, one file per reason | every exclusion is listed by ID, not just counted |
| `APM_locus/` | chr6 HLA/APM specifics | the locus the paper is about needs its own handling |
| `results/` | figures and summary tables used by the documentation | what the md files cite |
| `docs/` | subway-map diagrams of both workflows | drawn from `snakemake --rulegraph` |
| `earlier_tests/` | superseded analyses kept for provenance | shows what was tried before the current design |

## Files

| file | what it is | why it is here |
|---|---|---|
| [README.md](README.md) | entry point: what the pipeline produces, headline numbers with sources | start here |
| [HOW2RUN.md](HOW2RUN.md) | requirements and commands for both workflows, local and SLURM | neither runs from a fresh clone without setup |
| [INTRO.md](INTRO.md) | the one-page introduction, with citations | the paper's framing |
| [METHODS.md](METHODS.md) | every step: what, why, which script, which reference | the methods text, traceable to code |
| [LITERATURE.md](LITERATURE.md) | 26 papers mapped to the pipeline step each supports | justification per decision |
| [LEARNED_LESSONS.md](LEARNED_LESSONS.md) | why this design and not the obvious alternatives | isoforms, UTRs, kallisto vs STAR, degradation |
| [PHENOMENA.md](PHENOMENA.md) | the mappability phenomena tested, with numbers and links | the problem catalogue |
| [LIMITATIONS.md](LIMITATIONS.md) | what the pipeline does not do, or measures with a known bias | bounded, with a source per entry |
