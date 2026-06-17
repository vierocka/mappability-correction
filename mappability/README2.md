# Mappability correction — rationale and design

## The bias problem

At loci where multiple highly similar sequences exist in the genome — pseudogene-dense regions, tandem gene families, repetitive elements — RNA-seq quantification tools fail in opposite directions depending on which reference and tool are used.

**Signal loss — STAR + featureCounts, full genome reference.**
Reads from a functional gene may align equally well to nearby pseudogene loci. featureCounts discards multimapping reads by default, causing systematic underestimation. The problem is most severe in gene-body exons where pseudogene sequences are nearly identical to the functional allele; UTR and flanking regions carry more divergent SNPs and are less affected.

**Signal inflation — Kallisto / Salmon, MANE Select reference.**
MANE Select contains one representative transcript per gene and no pseudogenes. Reads from expressed pseudogene transcripts are forced onto the most k-mer-compatible target — the functional allele. The EM algorithm cannot redistribute them because pseudogene sequences are absent from the reference. The result is inflated counts.

**Deduplicated CDS — intermediate case.**
After removal of exact CDS duplicates, pseudogene sequences that differ from the functional allele by even a few SNPs are retained. Inflation is attenuated relative to the MANE case because more competing sequence is present, but signal loss is also less severe than with the full genome.

**Short transcripts — mappability floor.**
Transcripts too short to generate any uniquely mappable read have a uniqueness factor of zero regardless of reference configuration. This is a baseline effect independent of pseudogene competition.

**UTR smoothing in pseudoalignment.**
The MANE RNA reference includes full transcripts with 5′ and 3′ UTRs. UTR sequences diverge faster between paralogs than CDS and carry more unique k-mers. Reads pseudoaligning to UTR regions are assigned confidently even when CDS k-mers are shared across pseudogenes. As a result, Kallisto against MANE RNA overestimates uniqueness relative to Kallisto against the deduplicated CDS reference, which lacks UTRs. The difference between the two Kallisto estimates isolates the UTR smoothing effect quantitatively.

**Intron retention and splicing defects.**
In cancer, reads with intronic sequence arise from intron retention, spliceosome mutations (SF3B1, U2AF1, SRSF2 are among the most commonly mutated splicing factors in haematological malignancies and solid tumours), and transcription read-through. These reads fail to pseudoalign against spliced-only references, producing signal loss orthogonal to the pseudogene inflation problem. Both effects can co-occur at the APM locus.

---

## Affected gene categories

The mappability problem is not specific to the APM locus. Several categories of protein-coding genes are systematically mis-quantified by standard pipelines.

**HLA class I genes.** The APM locus on chromosome 6p21 is the most polymorphic region of the human genome. Classical HLA class I (HLA-A, -B, -C) are surrounded by pseudogenes and paralogous sequences across multiple alt loci. Both signal loss (full genome reference) and signal inflation (MANE reference) are well documented here.

**Ribosomal protein genes (RPL/RPS family).** Approximately 2,000 ribosomal protein pseudogenes are annotated in the human genome — among the highest pseudogene densities of any gene family. Exons are often too short for reads to map exclusively to the functional locus. Ribosomal protein genes are also commonly used as normalisation reference genes, making this a calibration problem with downstream consequences for expression comparisons.

**Ubiquitin genes (UBB, UBC, UBA52, RPS27A).** The core ubiquitin coding sequence is nearly identical across all four genes. The vast majority of reads mapping to the ubiquitin gene family are multimappers, rendering expression essentially unquantifiable by pipelines that discard non-unique alignments.

**KRAB zinc finger proteins (KZNF / ZNF family).** Several hundred annotated KZNF genes constitute the largest transcription factor family in the human genome. They evolved to recognise and silence specific classes of repetitive elements (LINE-1, SINE, ERV families) and share highly homologous zinc finger domain arrays. In cancer, endogenous retrovirus derepression drives upregulation of large numbers of KZNFs, but identifying which individual family members are affected is complicated by multimapping between paralogs.

**Olfactory receptor genes (OR family).** Approximately 400 functional OR genes and over 600 OR pseudogenes are annotated in the human genome, making this one of the pseudogene-densest gene families. Ectopic OR expression has been reported in several cancer types. Standard pipelines lose most of the signal in this family.

**Replication-dependent histone genes (HIST clusters).** Canonical histone genes are intronless — no splice junctions are available to help disambiguate reads between paralogs. Coding sequences are under strong purifying selection and are nearly identical across many copies within each cluster. Uniqueness factors for these genes are expected to be consistently low regardless of reference configuration.

**Protocadherins (PCDHA/B/G clusters).** Roughly 55 genes across three tandem clusters encode proteins with high within-cluster sequence similarity, particularly in the extracellular domains. Reported to be dysregulated in several cancers. Individual gene quantification within each cluster is unreliable.

**Cancer-testis antigens (MAGE, GAGE, SPANX, CT45, CT47 families).** Many cancer-testis antigen families reside in repetitive regions of the X chromosome and contain multiple members with high within-family sequence identity. Frequently activated in cancer and proposed as immunotherapy targets, but standard RNA-seq quantification of individual family members is severely affected by multimapping.

**Endogenous retroviruses and transposon-derived transcripts (HERV, LINE-1).** Thousands of HERV copies share sequence with each other and with expressed HERV-derived protein-coding genes. HERV-K (HML-2) is documented as expressed in melanoma, breast cancer, and germ cell tumours. LINE-1 elements are transcriptionally reactivated in many cancers. Note: most HERV loci are not represented in MANE Select; the uniqueness factors here cover only HERV-derived protein-coding genes present in the MANE annotation.

---

## Reference configurations

Three references are used for STAR, representing a hierarchy from most to least competing sequence context:

| Reference | File | Sequences | Content |
|---|---|---|---|
| Full genome | GCF_000001405.40_GRCh38.p14_genomic.fna | 705 | Primary chromosomes, alt loci, unplaced scaffolds |
| Deduplicated CDS | GCF_000001405.40_GRCh38.p14_cds_from_genomic.dedup.fna | 93,088 | CDS only (no UTRs); 53,249 exact duplicates removed with seqkit rmdup |
| MANE Select RNA | MANE.GRCh38.v1.5.ensembl_rna.fna | 19,437 | Full spliced transcripts including UTRs; one per gene |

The 73,651 sequences in the deduplicated CDS beyond MANE coverage represent CDS from non-MANE isoforms, pseudogenes, and other annotated genes with unique sequence. These are the competing sequences that drive multimapping at loci with complex annotation geometry such as the HLA locus. MANE transcript CDS sequences are a subset of the 93,088 — they survived deduplication as representative sequences.

For Kallisto, two indices are used: MANE RNA and deduplicated CDS (the full genome is excluded; Kallisto is designed for transcriptome-scale references).

**Translation note.** The deduplicated CDS file is directly translatable — sequences begin at ATG, frame 1 throughout (NCBI convention). The MANE RNA file requires UTR trimming before translation; CDS coordinates are available from the MANE GFF. Some NCBI CDS entries include the stop codon in the sequence — verification before batch translation is advisable. Pseudogene CDS entries may carry internal stop codons or frameshifts and should be filtered by `gene_biotype=protein_coding` before use.

---

## Analysis design

### STAR — uniqueness factor

Reads are simulated and aligned back to each reference with STAR. The uniqueness factor (UF) is the fraction of simulated reads that recover MAPQ = 255 (unambiguously mapped), computed per transcript from R1 reads only, independent of any downstream counting tool.

**v1–v6 (43 variants) — MANE transcript-based reads.** Reads simulated from all 19,437 MANE Select transcripts or from genomic exon ± 50 bp intronic flanks. UF is per MANE transcript. Parameter sweep covers reference, read length (75/100/150/200 bp), mode (PE/SE), and read source.

In v2 settings (dedup CDS reference), reads from all 19,437 MANE transcripts — including non-coding RNAs such as snRNAs — are simulated from the MANE RNA source, but the dedup CDS reference contains only protein-coding sequences (93,088 CDS entries, all [gbkey=CDS]). Reads originating from non-coding transcripts cannot align to this reference because no matching sequence exists there. These transcripts appear in the v2 output with UF = 0. This reflects reference scope — the CDS reference excludes non-coding genes entirely — not sequence-identity multimapping. Do not interpret v2 UF = 0 for snRNAs or lncRNAs as equivalent to the multimapping-driven UF = 0 observed in MANE-reference settings.

**v7 (10 variants) — deduplicated CDS reads.** Reads simulated directly from the 93,088 deduplicated CDS sequences (protein-coding only; all [gbkey=CDS]) with a step-1 sliding window. UF is per CDS sequence. Two alignment targets: full genome (v7_a–f, v7_i–j) and dedup CDS self-map (v7_g–h). Non-coding transcripts are not present in the CDS source FASTA — no reads are ever generated for them. The simulation scope is entirely protein-coding; UF = 0 entries for non-coding features in the downstream annotation are out-of-scope artefacts, not mappability measurements.

**v8 (3 variants) — proper PE simulation.** Reads simulated from MANE transcripts with realistic 300 bp insert geometry: R1 from the 5′ end of the fragment, R2 from the reverse complement of the 3′ end. This matches real Illumina PE 150 bp library structure. v8_a aligns to the full genome, v8_b to dedup CDS, v8_c to the MANE transcriptome. v8_b and v8_c reuse the reads generated by v8_a. The same protein-coding scope caveat as v2 applies to v8_b: non-coding transcript reads cannot align to the dedup CDS reference.

### featureCounts — counting-level cross-validation

featureCounts is applied to existing STAR BAM files — no new alignment needed. This directly tests the TCGA/GDC standard pipeline and shows how the MAPQ-based uniqueness factor translates into actual count-level distortion.

### Kallisto — pseudoaligner cross-validation

Three read sources tested against both indices:

| Read source | Expected outcome |
|---|---|
| MANE spliced RNA (v5 reads) | Inflation at paralogy-dense loci |
| Exon-flank (intronic flanks, v5 reads) | Failure to pseudoalign; signal loss |
| Dedup CDS (v7 reads) | CDS-level uniqueness in pseudoalignment |

Comparing MANE RNA vs deduplicated CDS index for the same read set isolates the UTR smoothing contribution to the uniqueness estimate. Script: `Map_scripts/kallisto_mappability_crossval.sh`.

---

## Expected bias by configuration

| Reference | Tool | Expected bias at HLA/APM locus |
|---|---|---|
| Full genome | STAR + featureCounts | Underestimation — multimappers discarded |
| Full genome | STAR (MAPQ=255 UF) | UF directly quantifies signal loss |
| Deduplicated CDS | STAR + featureCounts | Underestimation, less severe than full genome |
| Deduplicated CDS | Kallisto | Overestimation — pseudogene reads collapse onto functional allele; no UTR rescue |
| MANE RNA | STAR + featureCounts | Moderate; some multimapping due to splice variants not in MANE |
| MANE RNA | Kallisto | Overestimation — UTR rescue further inflates signal |

---

## Gene classification framework

Combining the STAR uniqueness factor (v1–v7) with the Kallisto inflation estimate (ratio of abundance between MANE RNA and dedup CDS indices) defines four operational gene classes:

| Class | STAR UF (full genome) | Kallisto signal (MANE vs dedup CDS) | Action |
|---|---|---|---|
| **Always OK** | ≈ 1.0 | No inflation | Any standard pipeline is reliable |
| **Apply correction** | < 1.0, > 0 | No or mild inflation | featureCounts underestimates; divide observed count by UF |
| **Prefer pseudoaligner** | ≈ 1.0 | Inflated | STAR + featureCounts correct; Kallisto/Salmon overestimate |
| **Always wrong** | ≈ 0.0 | Inflated | Both tools biased in opposite directions; unrecoverable from short reads |

Expected class assignments:
- *Always OK* — single-copy housekeeping genes in unique genomic regions
- *Apply correction* — partial pseudogene overlap; most HLA genes with dedup CDS reference; many KZNF genes
- *Prefer pseudoaligner* — genes where competing sequences are not transcribed and absent from MANE
- *Always wrong* — canonical histones, ubiquitin, ribosomal proteins with high pseudogene density, olfactory receptors, HLA with full genome reference

The **apply correction** class is the primary target of this project. The **always wrong** class marks genes where proteomics or long-read sequencing is the only reliable quantification route.

---

## Correction validity and limits

**Healthy tissue.** In a non-malignant cell, pseudogenes at the HLA/APM locus are largely transcriptionally silent. Multimapping occurs because functional-gene reads share sufficient sequence identity with pseudogene loci to receive sub-maximal MAPQ — not because pseudogenes produced competing reads. Under this assumption the correction is valid: true\_count ≈ observed\_count / UF.

**Allele mismatch at the population reference.** HLA is the most polymorphic locus in the human genome. Exons 2 and 3 of HLA class I (encoding the peptide-binding groove) differ substantially between alleles. If the sample carries an allele absent from GRCh38, reads from divergent exons fail to map or receive low MAPQ — not due to pseudogene competition but due to reference mismatch. The UF computed from the reference allele then underestimates the true uniqueness of the sample's actual allele, and the correction overcorrects.

**Personalised genome approaches** resolve the allele-mismatch component. Once the individual's actual HLA alleles are included in the reference, divergent exons map uniquely and the UF is higher. Personalised references are more accurate for individual-level analysis but require per-sample HLA typing and cannot be applied retroactively to archived cohorts. The UF approach is complementary: cohort-applicable without per-sample typing, at the cost of slight overcorrection for samples with alleles highly divergent from GRCh38.

**Tumour tissue.** The correction assumption breaks down in cancer: pseudogene reactivation (cancer-testis antigens, endogenous retroviruses), HLA loss of heterozygosity, allele-specific silencing, and somatic mutations all alter the expected read distribution. In these contexts the UF derived from healthy-cell assumptions is inappropriate. Consistent with RNA underestimation from multimapping, a proteogenomics study found HLA-I protein significantly overexpressed relative to RNA in 41–63% of cancer samples ([Gao et al., *J Proteome Res* 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10629274/)). Protein-level HLA measurement (immunopeptidomics) or long-read sequencing is required for reliable allele-specific quantification in tumour samples.

---

## Proteomics does not bypass the ambiguity

Recommending protein-level measurement as a fallback assumes that mass spectrometry resolves the ambiguity. It often does not, for reasons that mirror the RNA-seq biases:

- *Shared peptides / protein inference.* Most tryptic peptides from HLA class I, ribosomal proteins, and ubiquitin are shared between paralogs. A peptide from the conserved α3 domain of HLA-A is identical in HLA-B and HLA-C. The peptide-to-protein assignment has no unique solution for shared peptides — the direct analogue of RNA multimapping ([Nesvizhskii & Aebersold, *Mol Cell Proteomics* 2005](https://pubmed.ncbi.nlm.nih.gov/16009968/)).

- *Database dependency.* Peptide-spectrum matching requires the correct protein sequence to be present in the search database. If the sample's HLA allele is absent, allele-specific peptides in the variable α1/α2 domains go unmatched. The choice of search engine alone shifts allele-specific HLA peptide identification by 30–70% ([Parker et al., *Mol Cell Proteomics* 2021](https://pubmed.ncbi.nlm.nih.gov/34303857/)); allele-matched spectral libraries improve detection 2–3-fold ([Pak et al., *Mol Cell Proteomics* 2021](https://pubmed.ncbi.nlm.nih.gov/33845167/)). This is the protein-level analogue of allele-mismatch in RNA-seq.

- *FDR inflation from large databases.* Using an overly large search database (all isoforms, variants, non-canonical ORFs) increases chance matches and distorts the false discovery rate — the protein analogue of alignment specificity loss with an overly complex genomic reference.

- *Pseudogene-derived and non-canonical peptides.* At least 1,546 pseudogenes produce detectable protein products by mass spectrometry ([Vasylieva et al., *J Proteome Res* 2024](https://pubmed.ncbi.nlm.nih.gov/39486438/)). Non-canonical proteins from cryptic translation are absent from standard databases yet populate the immunopeptidome ([Ruiz Cuevas et al., *Cell Reports* 2021](https://pubmed.ncbi.nlm.nih.gov/33691075/)). These are either missed or misassigned to the canonical gene — the protein-level analogue of Kallisto inflation.

The conclusion is that proteomics shifts the ambiguity problem from k-mer space to peptide-spectrum space. A personalised approach — allele-typed reference database and sample-matched spectral library — is needed for reliable HLA protein quantification, just as a personalised genome is needed for RNA-seq.

---

## Implications for prediction tools

Any model using RNA-seq expression of APM/HLA genes as input inherits the quantification bias:

- **Immunotherapy response prediction** — HLA-I expression is a standard predictor of checkpoint inhibitor response. If STAR + featureCounts systematically underestimates HLA expression with a full genome reference, models trained on TCGA will miscalibrate when applied to cohorts quantified differently.
- **Immune deconvolution** (TIMER, CIBERSORT, xCell, MCP-counter) — deconvolution signatures for cytotoxic T cells and NK cells include HLA class I and APM components. Signature weights calibrated on one quantification platform will be off on another.
- **Antigen presentation scoring** — combined scores from TAPBP, TAP1/2, PSMB8/9, and HLA class I subunits classify tumours as antigen-presentation competent or deficient. Genes in the *always wrong* class will push this score in opposite directions depending on the tool.
- **Cross-cohort reproducibility** — TCGA/GDC uses STAR + RSEM + featureCounts; GTEx uses STAR + RSEM; some studies use Salmon or Kallisto. The same gene can appear up- or downregulated across cohorts purely as a tool artefact at loci with complex annotation geometry.

The gene classification produced here provides a principled filter: features in the *always wrong* class should be excluded from model training or replaced with protein-level measurements; features in the *apply correction* class should be corrected before training; only *always OK* features are directly comparable across pipelines without adjustment.

---

## Relation to prior work

The individual components have precedents; the specific combination appears to be novel.

**Read-back simulation for mappability assessment** is established: [BlackOPs (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3799449/) simulates reads from the reference, aligns them back, and flags mismapping positions. The approach here follows the same logic but produces a continuous per-transcript uniqueness factor rather than a positional blacklist, and targets RNA-seq rather than WGS.

**Genome vs. transcriptome mapping comparison** has been described ([PLOS One 2014](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0101374)), and systematic underestimation of HLA and immune gene expression by standard pipelines is documented ([Genome Biology 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0734-x); [Sci Rep 2023](https://www.nature.com/articles/s41598-023-41085-6)). These studies identify the problem but do not compute reference-specific correction factors.

**Multimapping correction tools** such as [CoCo](https://academic.oup.com/bioinformatics/article/35/23/5039/5505419) and [MGcount](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9158184/) operate downstream of alignment and do not assess how reference content itself shapes the mappability landscape.

**HLA-specific pipelines** ([HLApers](https://link.springer.com/protocol/10.1007/978-1-0716-0327-7_7), [HLAProphet](https://www.biorxiv.org/content/10.1101/2023.01.29.526142.full.pdf)) address quantification errors by building personalised references matched to the sample's HLA type — a complementary strategy that corrects the reference rather than characterising the distortion empirically.

What distinguishes the present approach: (i) MANE Select as the transcript set, (ii) three reference configurations compared within one framework, (iii) systematic parameter sweep over read length, PE/SE mode, and read source, (iv) the APM/HLA locus as primary target with the goal of diagnosing the direction and magnitude of bias per configuration, (v) Kallisto cross-validation to test the opposing inflation bias, (vi) featureCounts cross-validation to connect the MAPQ-based uniqueness factor to actual count-level output.
