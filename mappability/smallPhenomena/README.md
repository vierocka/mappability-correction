# smallPhenomena — Case Studies in RNA-seq Mappability

Each subfolder is a self-contained analysis that loads data directly from the
parent directories (`../results/`, `../kallisto_results/`, `../../Ref/`, etc.)
without depending on outputs from other subfolders.

Two shared libraries live in this folder and are imported by all scripts:

| Library | Purpose |
|---------|---------|
| `lib_check.py` | Read-length consistency checks before any statistics; documents that UF/TPM are ratio-based (cross-run comparison valid), while `est_counts` must never be compared across runs |
| `case_study_lib.py` | Shared loader and classifier for case studies 08A/08B: ENST→gene map from MANE FASTA, STAR/kallisto run metadata for all read lengths (75/100/150/200 bp, SE/PE, RNA/exflank, proper-pair), ambiguity classifier |

## Run order

All scripts are independent.  Run any in any order:

```bash
cd 01_isoform_redundancy     && python analyze.py
cd 02_overlapping_annotation && python analyze.py
cd 03_repetitive_pseudogene  && python analyze.py
cd 04_multicopy_coding       && python analyze.py
cd 05_local_duplications     && python analyze.py
cd 06_unique_genes           && python analyze.py
cd 07_sparse_celltypes       && python analyze.py
cd 08A_immune_evasion        && python analyze.py
cd 08B_HLA_MHC               && python analyze.py
```

`01_isoform_redundancy` runs MAFFT (available on PATH) and may take 10–30 min
for the 500-gene homology batch.

Case studies 08A and 08B load **all** available STAR and kallisto runs (all
read lengths and modes) via `case_study_lib.py`, so they are slower than the
single-run analyses.  Each `article/` subfolder contains the deposited PDF
and any supplementary files for the corresponding paper.

---

## Ambiguity source framework

These analyses collectively populate a classification of *why* a gene is
non-uniquely mappable.  Each row in the table below links to the primary
analysis subfolder that provides evidence for that ambiguity class.

| Ambiguity source          | Example                                    | Section       |
| ------------------------- | ------------------------------------------ | ------------- |
| Isoform ambiguity         | Alternative CDS transcripts of same gene   | 01            |
| Paralog ambiguity (local) | PCDH cluster, OR gene families             | 04, 05        |
| Paralog ambiguity (distal)| HLA-A vs HLA-H, RPL paralogs               | 04            |
| Genomic overlap           | Antisense lncRNA over protein-coding gene  | 02            |
| Nested gene               | Small gene inside intron of larger gene    | 02            |
| Pseudogene ambiguity      | RPL gene vs its processed pseudogenes      | 03            |
| Readthrough / retrogene   | GENCODE readthrough_gene, retrogene tags   | 02            |
| Multi-copy identical loci | RNU1, RNU6, histone clusters               | 04            |

Case studies 08A and 08B provide concrete published-paper examples spanning
multiple ambiguity classes simultaneously:

| Ambiguity source          | Example                                    | Section       |
| ------------------------- | ------------------------------------------ | ------------- |
| Allelic ambiguity         | HLA-A/B/C alleles, thousands of known variants | 08B       |
| Paralog ambiguity (HLA)   | HLA-A vs HLA-E vs HLA-F; HLA-DRA vs HLA-DPA1 | 08A, 08B  |
| Pseudogene ambiguity (HLA)| HLA-DRB6 — pseudogene in MANE, no UF data | 08B           |
| Mixed-class gene panel    | GBP1 (GBP family), IFITM2 (IFITM family), RHOG (Rho GTPase), HLA-B, HLA-E | 08A |

Two ambiguity sources are **not** covered by the current data:
- **Repeat ambiguity** (LINEs, SINEs, HERVs): these require repeat-masking
  analysis; not in the GENCODE annotation.
- **Allelic ambiguity** (HLA alleles, SNP-containing regions): requires
  population-level allele databases (e.g., IPD-IMGT/HLA) and phased genotypes.
  Sections 08A and 08B illustrate where this matters in practice, but the UF
  data here uses the single MANE canonical transcript per gene (no allelic
  decomposition).

---

## Data sources used

| Data file | Used by |
|-----------|---------|
| `../../Ref/GCF_000001405.40_GRCh38.p14_cds_from_genomic.dedup.fna` | 01 |
| `<SELF_NONSELF>/Reference/GENCODE/gencode.v49.annotation.gtf.gz` | 01 (pseudogene counts) |
| `<SELF_NONSELF>/Reference/MANE/*_enst.txt` (kinase/tf/cell_cycle/extracellular/blacklist) | 01 (functional groups) |
| `../../Ref/MANE.GRCh38.v1.5.ensembl_rna.fna` | 02, 03, 04, 05, 06 |
| `../../Ref/MANE.GRCh38.v1.5.ensembl_genomic.gff` | 02 |
| `../results/cds_uniqueness_factors_selfmap_PE.tsv` | 01 |
| `../results/id_conversion_table.tsv` | 02, 03 |
| `../kallisto_results/mappability_dedup_cds_reads_PE/kallisto_vs_dedup_cds/abundance.tsv` | 01 |
| `../kallisto_results/mappability_MANE_RNA_PE/kallisto_vs_MANE/abundance.tsv` | 02, 05, 06 |
| `../sc3UTRs/results/gene_overlap_annotation.tsv` (18 MB, gitignored; rebuild with `../sc3UTRs/annotate_gene_overlaps.py`; per-phenomenon counts are committed as `gene_overlap_phenomenon_counts.tsv`) | 02, 03, 04, 05, 06 |
| `../sc3UTRs/results/sc3prime_v5g_MANE_RNA_PE_100.tsv` | 02, 03, 04, 05, 06 |
| `../sc3UTRs/results/sc3prime_v2e_dedupcds_RNA_PE_100.tsv` | 03 |
| `../sc3UTRs/results/celltype_impact_summary.tsv` | 07 |
| `../sc3UTRs/results/marker_flags.tsv` | 07 |
| `../sc3UTRs/results/marker_3prime_impact.tsv` (78 MB, gitignored; rebuild with `../sc3UTRs/intersect_markers.py`) | 07 |
| `../sc3UTRs/results/sc3prime_v*.tsv` (all MANE runs) | 08A, 08B |
| `../kallisto_results/mappability_MANE_*/abundance.tsv` (all runs) | 08A, 08B |

---

## Key questions addressed per section

### 01 — Isoform redundancy
- How many genes in the deduplicated CDS have multiple isoforms?
- When reads are simulated from the dedup CDS and mapped back to itself,
  what fraction of positions map uniquely? (~6% overall)
- Which isoform pairs are >95% identical (MAFFT)?
- Do genes with more isoforms accumulate more pseudogenes? (GENCODE v49)
- Do multi-isoform genes cluster in specific functional classes (kinase, TF,
  extracellular, cell cycle)?
- Per isoform: how many unique 31-mers guide the signal vs shared k-mers that
  cause EM redistribution? Computed for CDS-only and full MANE transcript
  (UTR contribution) to expose reference-geometry batch effects.
- Exon-flank reads against both indices as a simplified total-RNA model: reads
  spanning exon-intron boundaries fail pseudoalignment, approximating what
  happens to intron-retaining or pre-mRNA reads in total RNA-seq studies.

### 02 — Overlapping annotation (MANE level)
- What fraction of MANE-level protein-coding genes overlap with another
  annotated element (antisense, nested, readthrough)?
- How does UF distribute across phenomena?  Which is the most disruptive?
- For overlapping groups, how does kallisto redistribute reads among members?

### 03 — Pseudogene ambiguity
- Which protein-coding genes are co-located with a pseudogene(s)?
- Does the 3' UTR help distinguish the functional copy?
  (`delta_3p250 > 0` → UTR adds uniqueness)
- How often does the UTR rescue a gene that has zero overall UF?

### 04 — Multicopy genes with high coding homology
- Which PC gene families have multiple copies in a sequence cluster?
- For these genes, does the UTR rescue uniqueness?
- Which families are worst (UTR also lost or no UTR rescue)?

### 05 — Local duplications
- For local tandem duplicates, how do STAR (discard) and kallisto (EM
  redistribute) differ in read assignment?
- Does STAR discard reads that kallisto confidently assigns (inflation risk)?

### 06 — Truly unique genes
- How many MANE protein-coding genes are genuinely unique?
- Do they have UF ≈ 1.0 as expected?
- Are there any "unique" genes with unexpectedly low UF (repeats in gene body)?
- Cross-phenomenon comparison: unique vs all other categories.

### 07 — Sparse scRNA-seq cell type risk
- Which cell types (≤200 markers) are most vulnerable to marker signal loss?
- With 20% detection rate, how many reliable markers remain after removing
  genes with uf_3p250 < 0.8?
- CRITICAL: < 2 expected reliable markers; HIGH: < 5; MEDIUM: < 10

---

## Case studies from published papers

### 08A — Immune evasion gene panel (tuberculosis biomarker)

**Paper:** Shi et al. (2022) "Role of GBP1 in innate immunity and potential as
a tuberculosis biomarker." *Scientific Reports* 12:11097.
doi:10.1038/s41598-022-15482-2

**Deposited files:** `08A_immune_evasion/article/s41598-022-15482-2.pdf`

**Gene panel (11 genes):**
GBP1, HLA-B, ELF4, HLA-E, IFITM2, TNFRSF14, CD274, AIM2, CFB, RHOG, HORMAD1

These genes were identified as hub genes predictive of tuberculosis (TB) using
WGCNA on the GSE83456 dataset (92 TB patients, 61 healthy controls).  GBP1
itself was the primary biomarker; the panel is used here as a real-data example
of a published gene list that mixes multiple ambiguity classes:

- **GBP1** — GBP family (GBP1–7), highly similar paralogs, potential inflation
- **HLA-B, HLA-E** — MHC class I, allelic and paralog ambiguity
- **IFITM2** — IFITM family (IFITM1/2/3), local tandem duplicates on chr11
- **RHOG** — Rho GTPase superfamily, large distal paralog family
- **TNFRSF14** — TNF receptor superfamily, large paralog family
- **CD274 (PD-L1), AIM2, CFB, ELF4, HORMAD1** — may be uniquely mappable;
  included as within-panel controls

**Key questions:**
- Which of the 11 hub genes have compromised UF (< 0.8) in 10x-relevant 3' windows?
- Do longer reads (150/200 bp) rescue the ambiguous genes?
- Does kallisto inflate TPM for GBP1/IFITM2 relative to the true signal?
- Which genes are reliable anchors for expression comparison?

**Outputs:** `gene_classification.tsv`, `star_all_runs.tsv`,
`kallisto_all_runs.tsv`, `wide_uf_table.tsv`

---

### 08B — HLA/MHC antigen presentation pathway (ICB immunotherapy)

**Paper:** Gong & Karchin (2023) "Clustering by antigen-presenting genes
reveals immune landscapes and predicts response to checkpoint immunotherapy."
*Scientific Reports* 13:950. doi:10.1038/s41598-023-28167-1

**Deposited files:**
- `08B_HLA_MHC/article/s41598-023-28167-1.pdf`
- `08B_HLA_MHC/article/41598_2023_28167_MOESM1_ESM.docx` (Supplementary)
- `08B_HLA_MHC/article/41598_2023_28167_MOESM2_ESM.docx` (Supplementary)

**Gene panel (32 genes):**

| Group | Genes |
|-------|-------|
| MHC class I | HLA-A, HLA-B, HLA-C, HLA-E, HLA-F, HLA-G |
| MHC class II | HLA-DRA, HLA-DRB1, HLA-DRB5, HLA-DRB6, HLA-DQA1, HLA-DQB1, HLA-DQA2, HLA-DQB2, HLA-DPA1, HLA-DPB1, HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB |
| Antigen presentation | B2M, NLRC5, CIITA, TAPBP, TAP1, TAP2, PSMB8, PSMB9, PSMB10, PSME1, PSME2, ERAP1 |

The paper developed a 23-gene APM (antigen-presenting machinery) expression
signature and used hierarchical clustering to stratify melanoma patients by
immune landscape and ICB response (Liu anti-PD1, Van Allen anti-CTLA4 cohorts).
This panel is the canonical real-data example for allelic and paralog ambiguity
in the MHC locus (chr6p21.3):

- **HLA class I/II** — thousands of known alleles; MANE carries only one
  canonical transcript per gene; all allelic diversity is invisible to
  standard mapping
- **MHC chr6p21.3 density** — most gene-dense, polymorphic, repeat-rich region
  in the genome; genomic overlap and local duplication effects compound allelic
  ambiguity
- **HLA-DRB6** — `transcribed_unprocessed_pseudogene` (GENCODE v49, no MANE transcript,
  no UF data); included because the paper listed it — see
  `HLA_DRB6_pseudogene_as_biomarker.txt` for the full discussion.  The
  `HLApseudo_MAFFT/classII_DRB` alignment shows that DRB6 (715 bp) contributes
  to the 66% shared-window rate at 100 bp for HLA-DRB1 — its "expression" in
  bulk RNA-seq is a mixture of genuine pseudogene transcription and mappability
  leakage from functional DRB genes
- **Non-HLA genes (B2M, PSME*, PSMB*, etc.)** — expected to be uniquely
  mappable; serve as within-panel controls

**Key questions:**
- What is the UF for each MHC gene across all read lengths?
  Do longer reads substantially rescue class II genes?
- For the chr6 MHC locus genes, how does UF correlate with chromosomal
  position and local gene density?
- Which antigen-presentation machinery genes outside chr6 are reliable?
- What is the quantitative risk that low UF on HLA genes inflates (kallisto)
  or suppresses (STAR) their apparent expression in scRNA-seq?

**Outputs:** `gene_classification.tsv`, `star_all_runs.tsv`,
`kallisto_all_runs.tsv`, `wide_uf_table.tsv`, `mhc_locus_overview.tsv`

The `mhc_locus_overview.tsv` lists all chr6 genes from the panel sorted by
chromosomal start position — useful for visualising the APM cluster density
against UF values.

---

## Simulation design — read sources and indices

Three read sources were simulated; each was quantified against two kallisto indices:

| Read source | What was simulated | Biological proxy |
|---|---|---|
| **MANE RNA** | Sliding window over full spliced transcript (5'UTR + exons + 3'UTR) | Poly-A RNA-seq |
| **Exon-flank** | Each exon ± 50 bp of flanking intronic sequence | Simplified total RNA-seq / intron retention |
| **Dedup CDS** | Sliding window over coding sequence only (no UTRs) | CDS-level annotation redundancy |

| Index | Sequences | UTRs | Isoforms |
|---|---|---|---|
| **MANE RNA** | 19,437 MANE Select spliced transcripts | Yes | 1 per gene |
| **Dedup CDS** | 93,088 RefSeq CDS entries after exact deduplication | No | All annotated isoforms |

Full genome index excluded from kallisto — not designed for genome-scale references.
STAR was run against the full genome (v1–v6) and against the dedup CDS (v7).

**Exon-flank + kallisto as a total-RNA proxy:**
Reads landing entirely within an exon pseudoalign normally; reads spanning the
exon-intron boundary have no match in the spliced index and fail pseudoalignment.
This approximates what happens when total RNA (with intronic content) is run through
a kallisto/Salmon pipeline with a spliced transcriptome reference.  The ±50 bp
flank is conservative — real total RNA has intronic reads proportional to intron
length — so measured effects are lower bounds on the real protocol mismatch.

---

## Glossary

| Term | Definition |
|---|---|
| **UF** | Uniqueness factor = fraction of simulated reads recovering MAPQ 255 (STAR unique alignment) |
| **uf_3p250** | UF computed only within the 3′-terminal 250 bp window (relevant for 3′-end scRNA-seq) |
| **TPM** | Transcripts Per Million — length- and library-size normalised; safe to compare across samples of the same read length |
| **est_counts** | Raw EM-estimated read counts from kallisto; NOT normalised for library size — never compare across runs |
| **PE / SE** | Paired-end / single-end sequencing mode |
| **L75 / L100 / L150 / L200** | Simulated read length in bp |
| **MANE** | Matched Annotation from NCBI and EMBL-EBI — one canonical transcript per gene |
| **dedup CDS** | Deduplicated RefSeq CDS FASTA — all protein-coding isoforms, exact duplicates removed |
| **selfmap** | Reads simulated from a FASTA and aligned back to the same FASTA; isolates annotation-internal ambiguity |
| **exon-flank** | Reads simulated from exon ± 50 bp intronic flanks in genomic coordinates |
| **k-mer (31-mer)** | Sequence substring of length 31; kallisto and Salmon default index k-mer size |
| **unique k-mer** | A k-mer present in exactly one transcript in the reference index; guides unambiguous read assignment |
| **entropy** | Normalised Shannon entropy [0–1] of TPM distribution across isoforms: 0 = all signal on one isoform, 1 = uniform spread |
| **EM** | Expectation-maximisation — kallisto/Salmon algorithm for redistributing multi-mapping reads |
| **MAPQ 255** | STAR flag for uniquely mapping reads (no equally scoring alternative alignment) |
| **prefix-match pseudogene** | Pseudogene assigned to a parent gene by gene_name prefix (e.g. RPL32P1 → RPL32); approximate, population-blind |
