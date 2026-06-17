# sc3UTRs — 3' End Mappability for Single-Cell RNA-seq

## Motivation

10x Genomics 3' scRNA-seq (Chromium v3/v3.1) captures only the last ~200–600 bp of each
transcript (the 3' UTR and often just the terminal exon). The overall uniqueness factor
(UF) computed by the main mappability pipeline covers the entire transcript length and
**does not reflect what 10x actually sequences**.

A gene can have high overall UF but critically low 3' UF, making it effectively
multi-mapping in 10x data. This is especially relevant for:
- Cell-type marker genes used in cluster annotation
- Gene families with conserved coding sequences but variable UTRs
- Genes that share 3' ends with pseudogenes

**The 3' end UF problem is known in the literature** (e.g., Melsted et al. 2019, Srivastava
et al. 2019, Bruning et al. 2022 on transcript-level quantification challenges) but is
rarely quantified per-gene for marker selection. This analysis provides that quantification.

---

## Main Questions

1. **Which genes have a shared or conserved 3' end?** For each simulated-read setting (read
   length 75–200 bp, single-end/paired-end, reference: full genome / MANE / dedup CDS), which
   genes show low UF when restricted to the last 250, 400, 500, 600 bp of their transcript?
   Low 3' UF is a consequence of 3' sequence sharing — between paralogs, recently duplicated
   family members, or pseudogenes that copied from the 3' end of their parental gene.

2. **Which cell-type markers are most affected?** Intersect the 3' UF results with CellMarker 2.0
   and PanglaoDB marker tables. Which cell types lose reliable markers in 10x 3' data?

3. **Do transcript isoforms of the same gene share their 3' ends?** Using all GENCODE v49
   isoforms (not just MANE Select), test whether the last W bp is more often identical between
   isoforms of the same gene than a randomly sampled internal window (controlling for general
   isoform similarity). Shared 3' ends mean 10x reads are gene-level unambiguous but lose
   isoform resolution.

4. **Do genes share 3' ends with pseudogenes?** A match between a protein-coding transcript's
   3' end and a pseudogene transcript means reads may be misassigned — not multimapping in the
   STAR sense, but lost to the pseudogene locus in alignment. Pseudogene 3' end sharing is a
   known source of expression underestimation.

5. **Can we form cross-gene 3' end clusters?** Standard scRNA-seq assigns reads to genes
   (sum-of-isoforms model). If multiple genes or gene families share a 3' sequence, they form
   a confusion cluster that undermines per-gene quantification regardless of how single reads
   map. Clustering by 3' end sequence reveals which genes are co-ambiguous.

---

## Reference consistency

| Reference | Version | Assembly | Chr naming | Used for |
|-----------|---------|----------|-----------|----------|
| GCF_000001405.40_GRCh38.p14_genomic.fna | NCBI RefSeq, ann. 2025-08 | GRCh38.p14 | NC_000001.11 | STAR genome index (v1–v8 scripts) |
| MANE.GRCh38.v1.5.ensembl_rna.fna | NCBI/Ensembl MANE v1.5 | GRCh38 | — (transcript sequences) | Transcript FASTA for v5, v8c |
| gencode.v49.transcripts.fa.gz | GENCODE v49 / Ensembl 115 | GRCh38.p14 | — (transcript sequences) | 3' end similarity analysis |
| gencode.v49.basic.annotation.gtf | GENCODE v49 | GRCh38.p14 | chr1, chr2... | Gene/biotype annotation for 3' analysis |

**GCF_000001405.40_GRCh38.p14 and GENCODE v49's GRCh38.p14 are the same assembly** — only
the chromosome naming convention differs (NC_000001.11 vs chr1). For the 3' end sequence
analysis we compare transcript sequences directly (no genomic coordinates), so naming
differences are irrelevant. Genomic positions CAN differ between RefSeq and GENCODE gene
models for the same assembly — this matters for coordinate-based analyses but not for the
sequence-comparison approach used here.

MANE v1.5 ENST IDs correspond to Ensembl ~109–111. GENCODE v49 = Ensembl 115. For the small
subset where transcript versions differ, the base ENST ID (without version) remains stable for
MANE Select transcripts. Cross-reference uses gene symbols (stable across both).

---

## Scripts

### `compute_3prime_uf.py`
Single-script per-BAM analysis.
- Input: STAR BAM + per-transcript UF TSV (from main pipeline) + read length
- For each transcript, counts reads (R1 only, MAPQ=255) originating from the last W bp
  (W = 250, 400, 500, 600)
- 3' window definition: read position `pos >= n_positions - 1 - W + 1`
- Denominator: `min(W, n_positions)` (whole transcript if shorter than W)
- Output columns: transcript_id, gene_symbol, transcript_len, read_len, label, n_positions,
  uf_overall, n_sim_3p{W}, n_unique_3p{W}, uf_3p{W}, delta_3p{W} (×4 windows)
- Sorted by delta_3p500 (most 3'-degraded transcripts first)

### `run_3prime_HPC.sh`
SLURM wrapper. Calls `compute_3prime_uf.py` for all v1–v8 settings that use MANE RNA
or exon-flank reads (~36 BAMs). Skips exon-flank BAMs for 3' analysis by design
(read name format `tid|ex_idx|pos` would give wrong transcript ID via rfind('|')
— only MANE RNA BAMs produce correct `ENST...|pos` read names).
Output: `results/sc3prime_{label}.tsv` per run.

### `compute_3end_similarity.py`
3' end sequence comparison across all GENCODE v49 isoforms.
- Input: `gencode.v49.transcripts.fa.gz` + `gencode.v49.basic.annotation.gtf`
- Analysis 1 (isoform_sharing): For each gene with ≥2 transcripts, computes the fraction of
  isoform pairs sharing identical last-W bp; compares to fraction sharing first-W bp (background).
  Tests whether 3' sharing exceeds what would be expected from general isoform similarity.
- Analysis 2 (pseudogene_xref): For each protein-coding transcript, checks whether its last-W bp
  exactly matches any pseudogene transcript. Flags the matching pseudogene gene.
- Analysis 3 (3end_clusters): Groups all transcripts by identical last-W bp sequence. Reports
  cluster size, gene/biotype composition, and flags clusters that mix protein_coding and
  pseudogene biotypes.
- Uses: only transcript sequences — no genomic coordinates, no alignment.
- Outputs: `3end_isoform_sharing.tsv`, `3end_pseudogene_xref.tsv`, `3end_clusters.tsv`,
  `3end_cluster_summary.tsv` (per window W, in `results/`)

### `intersect_markers.py`
Marker gene impact assessment.
- Input: `results/sc3prime_*.tsv` + CellMarker 2.0 XLSX + PanglaoDB TSV.GZ
- For each marker gene: looks up 3' UF across all settings
- Flags markers where `delta_3p500 < -0.1` (≥10% UF drop in the 500 bp window)
- Output: per-cell-type summary (how many of its markers are affected)
- Also outputs a per-gene table with marker cell type(s) and 3' UF metrics

### `download_markers.sh`
Documents download commands for marker tables.
Currently downloaded files:
- `CellMarker2_Human_cell_markers.xlsx` — CellMarker 2.0 (NAR 2023); human cell markers
- `PanglaoDB_markers_27Mar2020.tsv.gz` — PanglaoDB (2020 release)

---

## Comparison plan

| Comparison | Question | Script |
|-----------|---------|--------|
| uf_overall vs uf_3p{250,400,500,600} | Which genes share 3' sequences with other loci? Low 3' UF = conserved/shared 3' end | compute_3prime_uf.py |
| uf_3p500 vs read_len (75,100,150,200) | Longer reads extend into more unique flanking sequence — can they rescue shared 3' ends? | intersect_markers.py |
| genome vs MANE reference | Does mapping to full genome reveal more or fewer shared 3' loci than MANE self-map? | intersect_markers.py |
| isoform 3' sharing vs background | Is 3' convergence stronger than general exon sharing? | compute_3end_similarity.py |
| protein-coding vs pseudogene 3' match | Which genes have pseudogene 3' shadows? | compute_3end_similarity.py |
| 3' end clusters vs gene families | Do confusion clusters cross gene family boundaries? | compute_3end_similarity.py |
| affected markers per cell type | Which annotations are most at risk in 10x? | intersect_markers.py |
