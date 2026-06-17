# Kallisto pseudoalignment cross-validation

## Motivation

The main mappability pipeline (v1–v8) uses STAR to align simulated reads and counts how many map uniquely (uniqueness factor, UF). A natural question is whether a **k-mer-based pseudoaligner** (kallisto) behaves differently for bulk RNA-seq.

Signal loss in alignment-based quantification has multiple causes that a pseudoaligner may or may not address:

- **Gene multiplication at high sequence identity**: when a gene family has members with near-identical or identical coding sequences, reads are divided equally among all copies. This smooths per-gene expression and attenuates differential expression signals even when a true difference exists.
- **Pseudogene transcription**: if a processed pseudogene is transcribed and its sequence is identical or near-identical to the parental gene, reads from the parental locus may be attributed to the pseudogene and vice versa, distorting both expression estimates.
- **Reference scope**: does the choice of reference (full transcript with UTRs vs. coding sequence only) change which reads are considered unique, and by how much?

Three open questions drive this cross-validation:

1. **Can a k-mer pseudoaligner rescue signal that STAR loses to multimapping?**
2. **Does pseudoalignment inflate signal, and through which mechanism?**
3. **Which reference is more appropriate — full transcript (MANE) or coding-only (dedup CDS)?**

Note: the UTR-based disambiguation aspect is primarily relevant for 3' single-cell RNA-seq (see `sc3UTRs/`), where only the last 200–600 bp of each transcript is captured. For bulk RNA-seq the full read length is available and the coding sequence redundancy is the dominant concern.

---

## Simulated reads — origin and counts

Reads are **the same simulated libraries used in the main v1–v8 STAR pipeline** — not newly generated. Each simulated library places one read at every possible position along the source sequences (exhaustive tiling). The number of reads equals the number of tiling positions.

| Read source | Description | Read lengths | n reads (approx.) |
|---|---|---|---|
| **MANE_RNA** | One read per position of each MANE Select transcript (19,437 transcripts, ~3,600 bp avg) | 75, 100, 150, 200 bp | 68–71 M per library |
| **MANE_exonflank** | Reads spanning exon–intron boundaries (±flank); tests splice-junction sensitivity | 75, 100, 150 bp | 62–77 M per library |
| **dedup_cds_reads** | One read per position of each deduplicated CDS entry (93,088 entries, ~2,100 bp avg) | 75, 100, 150, 200 bp | 180–192 M per library |

SE libraries are fed to kallisto as `--single -l {read_len} -s 1` (fragment length = read length, sd = 1).  
PE libraries supply R1 + R2 as a pair; kallisto infers fragment length from insert size.

---

## Kallisto indices

| Index | File | n targets | Source |
|---|---|---|---|
| MANE | `kallisto_index_MANE.idx` | 19,437 | `MANE.GRCh38.v1.5.ensembl_rna.fna` (full transcript: 5'UTR + CDS + 3'UTR) |
| Dedup CDS | `kallisto_index_dedup_cds.idx` | 93,088 | `GCF_000001405.40_..._cds_from_genomic.dedup.fna` (coding only, deduplicated) |

Kallisto version: **0.50.1**, k=31 (default), no bootstraps.

---

## Scripts

**`../summarise_kallisto.py`**  
Reads all `run_info.json` and `abundance.tsv` files across every condition × index combination and produces two output tables. Parses CDS target IDs (`lcl|NC_..._NP_..._N`) via regex to map NP_ → ENST for cross-index comparison.

**`../summarise_kallisto_genes.py`**  
Reads the per-transcript inflation table and classifies each transcript by its direction and consistency across all 22 conditions. Flags known multi-copy gene families (olfactory receptors, RPL/RPS, HLA, histones, hemoglobin) and prints top inflated/deflated genes.

---

## Folder structure

```
kallisto_results/
├── mappability_MANE_RNA_{SE,PE}_L{75,100,150,200}bp/     # 8 dirs
├── mappability_MANE_exonflank_{SE,PE}_L{75,100,150}bp/   # 6 dirs
├── mappability_dedup_cds_reads_{SE,PE}_L{75,100,150,200}bp/  # 8 dirs
│
│   Each read-source directory contains:
│   ├── kallisto_vs_MANE/
│   │   ├── abundance.tsv     # per-ENST TPM and est_counts
│   │   ├── abundance.h5      # HDF5 (same data)
│   │   └── run_info.json     # alignment stats
│   └── kallisto_vs_dedup_cds/
│       ├── abundance.tsv     # per-CDS-entry TPM
│       ├── abundance.h5
│       └── run_info.json
│
└── kallisto_results.tar.gz   # original archive from cluster
```

22 read-source directories × 2 indices = **44 kallisto runs** total.

---

## Output files (in `../results/`)

| File | Rows | Description |
|---|---|---|
| `kallisto_run_summary.tsv` | 44 | Per-run: n_processed, p_pseudoaligned, p_unique (from run_info.json) |
| `kallisto_inflation.tsv` | 425,968 | Per-transcript × condition: tpm_mane, tpm_cds, ratio, direction |
| `kallisto_run_overview.tsv` | 22 | Compact per-condition table with direction counts |
| `kallisto_genes_summary.tsv` | 19,436 | Per-transcript classification across all conditions |

---

## What differs between the two summary scripts

| | `summarise_kallisto.py` | `summarise_kallisto_genes.py` |
|---|---|---|
| Input | Raw kallisto output (run_info.json + abundance.tsv) | Pre-built kallisto_inflation.tsv |
| Granularity | Per-transcript × condition | Per-transcript (aggregated across all 22 conditions) |
| Key output | `kallisto_inflation.tsv` — the raw comparison table | `kallisto_genes_summary.tsv` — classification + gene families |
| ID handling | Parses NP_ from CDS target IDs, maps to ENST | Reads gene symbols directly from inflation table |
| Run order | Must run first | Depends on inflation table existing |

---

## Results summary

### Pseudoalignment rates

| Read source | → MANE index | | → CDS index | |
|---|---|---|---|---|
| | %pseudoaligned | %unique | %pseudoaligned | %unique |
| MANE_RNA | **100%** | 95–97% | 49–53% | 13–16% |
| exon-flank | 96–100% | 92–96% | 41–49% | 13–14% |
| CDS_reads | 98–99% | 94–96% | **~100%** | **5–6%** |

The striking numbers:
- **MANE RNA → CDS**: only ~50% of reads pseudoalign. The other 50% originate from UTR regions absent in the CDS reference — they have no compatible k-mers.
- **CDS reads → CDS self-map**: 100% pseudoalign but only 5–6% are unique. This is the true redundancy of coding sequences — without UTR disambiguation, nearly all CDS reads are ambiguous across paralogs.
- **CDS reads → MANE**: 95–96% unique. UTRs provide the k-mers that disambiguate what the CDS cannot.

### Per-gene classification (across 22 conditions)

| Class | n transcripts | Interpretation |
|---|---|---|
| `no_CDS_signal` | 3,740 | Non-coding transcripts (NR_ loci) — absent from CDS index, TPM = 0 |
| `mostly_no_CDS_signal` | 227 | Predominantly absent; CDS entry may exist for a subset of conditions |
| `MANE_strongly_inflated` | 3,488 | MANE TPM > 2× CDS TPM in ≥80% of conditions |
| `MANE_inflated` | 10,065 | MANE TPM > CDS TPM in ≥60% of conditions |
| `CDS_strongly_inflated` | 119 | CDS TPM > 2× MANE TPM in ≥80% of conditions |
| `CDS_inflated` | 864 | CDS TPM > MANE TPM in ≥60% of conditions |
| `mixed` | 700 | Direction varies by read length, SE/PE, or read source |
| `consistent` | 233 | Both indices agree; neither inflates |

**71% of transcripts are MANE-inflated** — UTR k-mers systematically rescue reads that are ambiguous in the coding space.  
**5% are CDS-inflated** — these genes have unique CDS k-mers but their MANE transcripts share UTR k-mers with other loci.

### Gene families

| Family | n | MANE-inflated | CDS-inflated | Note |
|---|---|---|---|---|
| Olfactory receptors (OR*) | 413 | 368 | 40 | Largest family; highly similar coding + UTR |
| Ribosomal protein L (RPL) | 53 | 25 | 23 | CDS redundancy strong; UTR partly rescues |
| Ribosomal protein S (RPS) | 46 | 30 | 11 | Same pattern as RPL |
| HLA | 21 | 20 | 1 | UTRs differ enough to rescue most |
| Histones H2A/H2B/H3/H4 | 84 | 51 | 26 | Split; shortest transcripts, minimal UTR |
| Hemoglobin alpha (HBA) | 2 | 0 | 2 | HBA1/HBA2: unique CDS, shared UTR |
| Tubulin α/β | 20 | 20 | 0 | UTR fully rescues |
| Calmodulin (CALM) | 7 | 2 | 3 | CDS nearly identical; partial UTR |

### Interpretation

The k-mer pseudoaligner **does not solve** the mappability problem — it shifts it. STAR loses reads to multi-mapping; kallisto redistributes them using UTR k-mers, inflating genes whose UTRs happen to be unique even when their CDS is shared. The CDS-only index reveals the underlying coding redundancy more honestly (~5–6% unique) but loses all non-CDS signal. For expression quantification of multi-copy families, both approaches have systematic biases that must be disclosed.
