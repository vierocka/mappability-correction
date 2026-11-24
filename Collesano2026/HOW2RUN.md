# How to run

Two Snakemake workflows. `server` turns FASTQ into corrected transcript TPM; `local` turns that
into the 9-mer source repertoire. They are separate because only the first needs a cluster.

## What each one needs before it will run

| | `server/Snakefile` | `local/Snakefile` |
|---|---|---|
| input data | paired-end FASTQ per run (fetched from SRA by the `fetch` rule, or local paths in `config/samples.tsv`) | the `server` pipeline's output under `results_server/` |
| references | `server/helpers/download_references.sh refs` — MANE v1.5 (RNA, genomic GFF, protein), Ensembl GRCh38 primary assembly, IPD-IMGT/HLA 3.63.1, TruSeq3 adapters | none |
| memory | ~40 GB for the STAR index of the primary assembly (the `funique` arm only; see below) | ~8 GB |
| other | conda (`--use-conda`), ~16 threads | — |

Neither will run on a laptop from a fresh clone: `server` needs the references downloaded and the
full FASTQs present, and `local` needs `server` to have finished first. A dry run (`-n`) checks
only that the DAG connects and that inputs exist — it does **not** prove the rules work.

**The `funique` arm normally does not run.** `resources/funique/` ships the precomputed uniqueness
tables for 75 bp and 100 bp reads, so Snakemake skips `funique_genome_index`/`funique_simulate`
unless they are missing. Only delete them if you intend to recompute, and then the 40 GB STAR
index is required.

## Local

```bash
cd Collesano2026
snakemake -s server/Snakefile --configfile config/server.yaml --use-conda --cores 16 -n   # dry run
snakemake -s server/Snakefile --configfile config/server.yaml --use-conda --cores 16
snakemake -s local/Snakefile  --configfile config/local.yaml  --use-conda --cores 4
```

## Running on an HPC cluster

```bash
cd Collesano2026
bash server/helpers/download_references.sh refs
snakemake -s server/Snakefile --configfile config/server.yaml --use-conda -n            # dry run
snakemake -s server/Snakefile --configfile config/server.yaml --use-conda \
  --executor slurm --jobs 20 --default-resources slurm_account=<account> mem_mb=16000   # Snakemake 8 slurm plugin
# results back: rsync -av <hpc>:<path>/Collesano2026/results_server/summary ./results_server/
```
