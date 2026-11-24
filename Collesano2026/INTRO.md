# Introduction

MHC class I molecules present a sample of the intracellular proteome to CD8⁺ T cells, and this self-immunopeptidome
defines what T cells learn to tolerate. Presentation is dynamic and strikingly inefficient: less than 0.1 % of the
peptides that survive intracellular degradation are presented (Yewdell 2003). Much of the supply comes not from
mature, long-lived proteins but from defective ribosomal products: more than 30 % of newly synthesised proteins
are degraded by proteasomes within minutes of synthesis (Yewdell 1996; Schubert 2000). Peptides also arise during
the pioneer round of translation of mRNAs destined for nonsense-mediated decay, and from pre-spliced RNAs in the
nucleus (Apcher 2011, 2013). Blocking cap-dependent translation suppresses full-length protein synthesis but only
partly affects antigenic peptide production (Apcher 2011). Such products are made, sampled and degraded within
minutes and never accumulate in the steady-state proteome. Reference-based shotgun proteomics therefore
under-represents exactly this source. It also cannot assign peptides shared by multicopy genes and paralogs to a
single protein (Nesvizhskii & Aebersold 2005). Bulk RNA-seq measures the transcripts that feed this supply,
including short-lived, decay-directed mRNAs. That makes transcript abundance, rather than protein abundance, the
appropriate proxy for antigen supply.

Proteogenomic studies support transcript abundance as a necessary but not sufficient determinant. Class I
peptides derive preferentially from actively transcribed, evolutionarily conserved genes. Yet 41 % of expressed
protein-coding genes yield no detectable peptides, and the HLA allotypes contribute unevenly
(Pearson 2016; Kubiniok 2022). Proteomics offers no simple alternative: the moderate mRNA–protein correlation
partly reflects protein measurement error (Upadhya 2022). Active NMD restricts the presentation of non-canonical,
frameshift-derived peptides, and inhibiting it reshapes the cancer immunopeptidome (Vendramin 2026). Which
transcripts are counted, and how, thus shapes any expression-weighted model of self.

Short-read quantification is least reliable precisely where it matters for antigen presentation. Mappability
varies widely between gene classes, and discarding multimapped reads underestimates paralog expression
(Derrien 2012). Twelve common methods misestimate hundreds of disease-relevant genes in multicopy families
(Robert & Watson 2015). Each strategy for multimapped reads trades bias for dilution (Deschamps-Francoeur 2020).
Pseudogenes compound the problem: many are transcribed and some are functional or translated
(Cheetham 2020; Lin 2025), yet most share long stretches of sequence with their parent genes. The HLA region adds
extreme polymorphism. Reads from non-reference alleles align poorly, and personalised references recover them
(Aguiar 2019; Kang 2023). HLA genotypes can be inferred directly from RNA-seq (Boegel 2012; Szolek 2014;
Orenbuch 2020) using curated allele sequences (Barker 2023), and cell-line records can be cross-checked in
Cellosaurus (Bairoch 2018).

Here we build a deliberately simple self model. Each protein-coding gene is represented by its MANE Select
transcript, a single representative, well-supported transcript (Morales 2022). HLA-A, -B and -C are replaced by
the sample's own alleles. Reads are aligned to this personal transcriptome (Dobin 2013; Li 2009) and multimapped
reads are kept as fractions. Expression is corrected for sequence uniqueness only where a pseudogene explains the
ambiguity. Every unique 9-mer then receives a net-sum expression: the sum of the TPM of all transcripts encoding
it. Applied to public RNA-seq (Leinonen 2011) of three cell lines with matched class I immunopeptidomes, the model
recovers 97.7–98.1 % of the eluted 9-mers. It thereby provides a reproducible, expression-annotated self
repertoire for testing expression as a variable in HLA-binding prediction.

*Full references: [LITERATURE.md](LITERATURE.md).*
