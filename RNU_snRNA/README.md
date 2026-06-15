# RNU snRNA — Mappability Analysis

## What this folder contains

| File | Description |
|---|---|
| `RNU_MANE.fna` | 23 RNU sequences extracted from MANE Select v1.5 RNA FASTA |
| `RNU_MANE_unique.fna` | Same sequences, identical copies collapsed (seqkit rmdup → 17 sequences) |
| `RNU1_family.fna` … `RNU6_family.fna` | Per-family FASTAs for alignment (4 multi-member families) |
| `RNU*_family.mafft.fna` | Per-family MAFFT alignments (`--globalpair --maxiterate 1000`) |
| `RNU_UF_summary.tsv` | Uniqueness factors for all 23 RNU genes across 4 MANE-reference simulation settings |
| `RNU_homology.pdf/png` | Main figure: identity heatmaps (RNU1 and RNU6 families) + UF bar chart |
| `extract_RNU_sequences.py` | Extraction and per-family FASTA generation |
| `run_mafft.sh` | Per-family MAFFT execution |
| `plot_RNU_homology.py` | Figure generation |
| `RNU_UF_summary.py` | UF summary table generation |

---

## Background: what snRNAs do

Small nuclear RNAs (snRNAs) are short, non-coding RNAs of 63–191 bp that form the catalytic core of the spliceosome. The major spliceosome contains U1, U2, U4, U5, and U6 snRNAs; the minor (ATAC) spliceosome contains U4ATAC, U6ATAC, U11, and U12. They are encoded at multiple genomic loci, do not have introns, and are not polyadenylated. Their accession prefix in RefSeq is `NR_` (non-coding RNA), unlike protein-coding transcripts which use `NM_`.

---

## Why MANE Select contains multiple copies

MANE Select (Matched Annotation from the NCBI and EBI) includes one representative transcript per gene per locus. snRNA gene families have expanded via retrotransposition, resulting in dozens to hundreds of genomic loci with nearly identical sequences. In GRCh38, there are >180 U1 pseudogene/gene loci and >965 U6 loci genome-wide.

MANE Select is correct to include all four RNU1 and five RNU6 loci that have curated gene status — they are distinct chromosomal loci with separate regulatory contexts. The redundancy within MANE is biologically accurate.

---

## The zero-signal problem

### Simulation setup

Using the `mappability-correction` pipeline, we simulated reads from the MANE Select RNA sequences and mapped them back to the MANE Select reference with STAR:

| Setting | Read type | Length | Reference |
|---|---|---|---|
| v5_a | SE MANE RNA | 75 bp | MANE Select |
| v5_a | SE MANE RNA | 100 bp | MANE Select |
| v5_d | SE MANE RNA | 200 bp | MANE Select |
| v5_b | PE exon-flank | 100 bp | MANE Select |

### Result

**RNU1-1, RNU1-2, RNU1-3, RNU1-4**: all four copies are 100% identical (164 bp). Every simulated read maps equally to all four loci → uniqueness factor (UF) = 0.000 in all settings. In the SE 200 bp setting, all are too short (164 < 200 bp) to simulate any reads at all.

**RNU6-1, RNU6-2, RNU6-7, RNU6-8, RNU6-9**: five near-identical copies (107 bp, 99.1–100% pairwise identity). Same result — UF = 0.000. All are also too short for SE 200 bp.

**RNU6ATAC**: the ATAC spliceosome U6 variant (126 bp). Despite carrying "RNU6" in its name, it is sufficiently diverged from the standard U6 copies (~70% pairwise identity) — UF = 1.000. It maps uniquely.

**RNU7-1**: 63 bp. Too short to simulate any reads with SE 75/100/200 bp settings. With PE exon-flank 100 bp, reads are generated but map to hundreds of genome-wide U7 pseudogenes → UF = 0.000.

The pairwise identity heatmaps and UF bar chart are shown in `RNU_homology.pdf`.

---

## Consequences for RNA-seq quantification

### Standard featureCounts / htseq-count

Reads simulated from these snRNAs multimap to all copies in the MANE reference. By default, featureCounts discards multimappers (`-M` not set). Result: **zero counts for RNU1-1/2/3/4 and RNU6-1/2/7/8/9**, regardless of true expression.

If these genes are actually differentially expressed, the difference will not be detected.

### Fractional counting (`featureCounts -O --fraction` or EM-based)

If fractional counting is enabled, each multimapping read contributes 1/N weight to each of the N targets. Since reads map equally to all N identical copies:

- Each locus receives 1/N of the signal
- Counts are artificially reduced by a factor of N
- All N copies receive identical counts (perfect artificial co-expression)
- Fold-change between conditions is numerically preserved, but detection power is reduced and no locus-specific information can be recovered

For a family with 4 identical copies (RNU1), each locus receives 25% of the true signal. For 5 copies (RNU6), each receives 20%.

### Kallisto / Salmon

Transcript-level quantifiers with EM handle multimapping probabilistically. For perfectly identical sequences, the EM algorithm cannot distinguish loci and distributes counts equally — same outcome as fractional counting. The effective signal per locus is 1/N of true expression.

---

## What would help

**Full genome reference instead of transcriptome**: When reads are simulated from the full genome including intergenic flanks (v1/v3 settings), the unique chromosomal context surrounding each locus can be captured. For RNU6 genes (107 bp), enough flanking sequence is unique per locus to recover partial signal (UF ≈ 0.88–0.91 in the full genome setting). For RNU1 genes (164 bp), the genomic flanks also overlap with additional pseudogene loci across the genome — UF remains near 0 even with full-genome mapping.

**Total RNA-seq**: Unlike poly-A selected RNA-seq, total RNA captures nascent pre-mRNA including introns. Intronic reads carry unique genomic flanking sequence (analogous to the exon-flank simulation). This is the most realistic path to locus-specific quantification of snRNAs from sequencing data.

**Long-read sequencing** (ONT, PacBio): With reads of 500+ bp, the unique flanking sequence per locus can be spanned even from poly-A depleted total RNA. Not yet standard for snRNA studies.

---

## Literature context

snRNA expression changes have been reported in several contexts:

**Alternative splicing regulation**: U1 and U2 snRNA levels influence splice site selection and global AS patterns. Dysregulation has been linked to neurodegeneration and cancer (reviewed in Matera & Wang, 2014, *Nat Rev Mol Cell Biol* 15:108–121).

**Cancer**: Recurrent U1 snRNA mutations (position 3, g.3A>C or g.3A>T) cause widespread cryptic splicing by altering 5' splice site recognition, detected in medulloblastoma, CLL, and hepatocellular carcinoma (Suzuki et al. 2019, *Nature* 573:142–146; Oh et al. 2020, *Nature* 566:173–179). These are somatic point mutations, not expression changes — RNA-seq detects them as mismatches, not as quantification signal.

**Polyadenylation and 3' processing**: U6 snRNA 3' end formation is Pol III-dependent; U1–U5 are Pol II transcripts with 3'-trimming. Neither class is polyadenylated → standard poly-A RNA-seq protocols systematically under-detect them.

**snRNA-seq**: Studies using Illumina sequencing to measure snRNA expression typically use total RNA or small RNA-seq protocols, and often collapse all copies of a family into a single pseudogene-aware reference, or report family-level rather than locus-specific counts. Locus-specific quantification of multi-copy snRNA families remains technically difficult.

---

## Key note on `RNU_MANE_unique.fna`

This file is the output of `seqkit rmdup` on `RNU_MANE.fna`. It contains the same MANE RNA sequences with identical copies collapsed to one representative — 17 sequences (6 removed: RNU1-2/3/4 and RNU6-2/8/9). It is **not** a coding-sequence extract; snRNAs have no CDS or UTR. Both files represent the complete, mature functional RNA molecules as defined by MANE Select.

---

## Summary

RNU1 (U1 snRNA) and RNU6 (U6 snRNA) genes in MANE Select are present as multiple near-identical or 100%-identical copies because they occupy distinct chromosomal loci. This biological reality makes locus-specific RNA-seq quantification impossible with standard poly-A RNA-seq and short-read transcriptome mapping. Simulated reads from these genes are indistinguishable by the aligner regardless of read length — uniqueness factor is 0 across all tested settings. Differential expression of individual RNU1 or RNU6 loci cannot be detected with conventional RNA-seq pipelines.
