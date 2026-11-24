# Mapping, counting and reference choice — consolidated notes

What the simulations in this project establish, in one place. Numbers are from
[`results/results_summary.tsv`](results/README.md), [`kallisto_results/`](kallisto_results/README.md),
[`smallPhenomena/`](smallPhenomena/README.md) and [`sc3UTRs/`](sc3UTRs/README.md).

## 1. The two failure directions

| | mechanism | direction |
|---|---|---|
| STAR + featureCounts, unique only | reads matching a pseudogene or paralog get MAPQ < 255 and are discarded | **loss** |
| kallisto / salmon on a reference without the competitor | EM cannot redistribute to a sequence that is absent, so reads collapse onto the functional copy | **inflation** |

Both act hardest on the same loci. A pseudoaligner does not solve the problem; it moves it from
multi-mapping to k-mer attribution.

## 2. Reference content decides the answer

Same reads, three references. Reads are simulated from the 19,437 MANE transcripts, single-end,
at 75/150/200 bp (`transcript_uniqueness_factors_{genomic,dedup_cds,MANE}_RNA_SE_L*.tsv`); each
cell is the range over those three read lengths. UF is per MANE transcript; "sequences" is the
size of the reference it was mapped against.

| reference | sequences | mean UF | UF = 1 | UF = 0 |
|---|---|---|---|---|
| full genome (GRCh38.p14, incl. alt/patches) | 705 | 0.89–0.90 | 11.0–14.2 k | 1,146–1,169 |
| dedup CDS (no UTRs) | 93,088 | 0.19–0.21 | 193–377 | 5.0–5.8 k |
| MANE transcriptome (one transcript per gene, UTRs kept) | 19,437 | 0.95–0.96 | 15.5–17.4 k | 141–173 |

The ordering is the whole story: **uniqueness is a property of the reference, not of the gene.** The same
transcript is "unique" against MANE and "ambiguous" against the genome, because the genome contains the
pseudogenes and alt contigs that MANE omits.

Caveat on the CDS column: reads simulated from all MANE transcripts include non-coding ones, which have no
target in a CDS-only reference and score UF = 0 for reasons of scope, not sequence. Do not read those zeros as
multi-mapping.

## 3. Coding sequence alone is not identifiable

Reads simulated from the deduplicated CDS set, aligned back:

| target | mean UF | global UF | UF = 0 |
|---|---|---|---|
| genome | 0.879–0.895 | 0.891–0.903 | ~5,900 |
| the CDS set itself | 0.109–0.110 | **0.056–0.057** | ~51,700 (56 %) |

Against itself the CDS reference resolves ~5.6 % of positions: more than half of all CDS entries share their
entire coding sequence with another entry. Genomic context — introns, flanks, UTRs — is what makes coding genes
distinguishable. This is the quantitative basis for keeping UTRs in any transcriptome reference.

## 4. UTRs are the discriminating sequence, and they also mislead

Simulated reads, kallisto, MANE (with UTRs) vs dedup-CDS index:

- CDS reads → CDS index: 5–6 % unique. CDS reads → MANE: 95–96 % unique.
- Unique k-mers: 3 % → 94 %.
- **71 % of transcripts are MANE-inflated** (MANE TPM > CDS TPM in ≥ 60 % of conditions); 5 % are CDS-inflated.

The 5 % matter: those genes have unique coding sequence but share UTR k-mers (HBA1/HBA2 is the clean example —
distinguishable CDS, shared UTR). So UTR rescue is not uniformly good; it is a second attribution channel with
its own failure mode.

1,922 genes are not rescued by UTRs at all: large families whose UTRs drift as slowly as their CDS (OR, ZNF,
KRTAP, PCDH, KIR ≈ 1,600), en-bloc duplications (PRR20, DEFB, AMY ≈ 200), replication-dependent histones with no
poly(A) 3′UTR (≈ 40) and Pol III ncRNA (RNU6 ≈ 17).

## 5. Read length and pairing

Over 75 → 200 bp against the genome, mean UF moves 0.886 → 0.903, and UF = 1 moves 11.0 k → 14.2 k. Longer reads
help, modestly and monotonically. Paired-end adds little over single-end at the same length; the gain comes from
total length, not from pairing.

Two artefacts to know: with PE reads against the MANE transcriptome the UF = 1 count collapses to ~0 while
UF ≥ 0.90 stays at 17.5 k, because fragment ends behave differently from read starts; and exon-flank reads
(simulating intron retention or pre-mRNA) against a spliced transcriptome drop to UF 0.68–0.74, since reads
crossing an exon boundary have no spliced target. Total-RNA libraries carry that penalty; poly(A) libraries
largely do not.

## 6. Why a gene is ambiguous — the source taxonomy

`sc3UTRs/annotate_gene_overlaps.py` labels every gene with one primary phenomenon, from GENCODE geometry plus
3′-end sequence clusters (transcripts sharing an **identical** terminal 250 bp; exact match, md5-keyed — note
the docstring in that script says "≥ 90 % identity", which the implementation does not do):

| phenomenon | n (GENCODE v49, all biotypes) |
|---|---|
| nested | 32,538 |
| genomic_overlap_antisense | 23,795 |
| unique | 17,627 |
| genomic_overlap_sense | 2,235 |
| local_tandem_dup | 1,286 |
| readthrough | 664 |
| retrogene | 391 |
| distant_paralog | 131 |
| pseudogene_shadow | 24 |

Only a minority of genes are geometrically clean. "Unique" is 22 % of annotated genes.

## 7. The action classes

Uniqueness (STAR, genome) and inflation (kallisto, MANE vs CDS) are independent axes, and together they say what
to do with a gene:

| class | UF | inflation | action |
|---|---|---|---|
| always OK | ≈ 1 | none | any pipeline is fine |
| apply correction | 0 < UF < 1 | none or mild | multiply the count by UF (never divide) |
| prefer pseudoaligner | ≈ 1 | inflated | STAR + featureCounts is right; kallisto overestimates |
| always wrong | ≈ 0 | inflated | both biased, opposite directions; needs long reads or protein |

Expected membership: housekeeping genes in unique regions; partial pseudogene overlap and most KZNFs; genes
whose competitors are untranscribed and absent from MANE; and histones, ubiquitin, high-pseudogene ribosomal
proteins and olfactory receptors respectively.

## 8. Correction validity

`corrected = observed × UF`, never divided. UF is the share of a locus's signal that is demonstrably its own,
so multiplying removes the share that cannot be attributed. The error is bounded by 1 − UF in either direction,
and its **sign depends on whether the competitor is transcribed**, which UF does not measure. In healthy tissue,
where MHC-region pseudogenes are largely silent, the correction removes real signal; in cancer, where they are
derepressed, it removes genuine inflation.

This is also why the correction must not fire when the competitor is *inside* the reference: those reads are
already split between records, and scaling again corrects twice.

## 9. Reference-source dependence

The `_nf` arm repeats the sweep against Ensembl GRCh38.116 primary assembly (nf-core/rnaseq default) instead of
NCBI RefSeq GRCh38.p14, to separate reference-choice effects from biology. Three distinct CDS references are
involved — genome-wide dedup CDS, MANE-restricted dedup CDS, and the genomic reference — and they are not
interchangeable; see `Map_scripts/README3.md` and `extract_dedup.log`.

Practical consequence: a uniqueness factor is only valid for the reference, annotation release, read length and
pairing it was computed on. Cross-cohort comparisons (TCGA with STAR+RSEM, GTEx with STAR+RSEM, others with
salmon or kallisto) inherit this, and a gene can appear up- or down-regulated purely as a tool artefact.

## 10. What is not covered

Repeat ambiguity (LINE/SINE/HERV) needs repeat-masking rather than annotation geometry. Allelic ambiguity needs
population allele databases and phased genotypes; everything here uses one canonical transcript per gene, so HLA
allelic diversity is handled separately by the personalised-reference route.
