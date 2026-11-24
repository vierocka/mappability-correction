# Literature behind the pipeline: claims and consistency check

Each paper's finding was checked against the pipeline step it supports. The PDFs are kept locally in `articles/`
(`SurnameYear_Journal.pdf`) and are **not** part of this repository (authors' and journals' rights).

Verdict key: **supports** = the pipeline follows the finding; **consistent, scoped** = no contradiction, but the finding
marks a boundary of the naive self model (stated in [LESSONS](LEARNED_LESSONS.md) and the [README limitations](../README.md#limitations)).

## Antigen presentation: why expression, and which expression

| paper | finding | pipeline step | verdict |
|---|---|---|---|
| Yewdell 1996 [Y96]; Schubert 2000 [S00] | DRiPs: > 30 % of newly synthesised proteins are degraded by proteasomes within minutes and feed MHC class I | expression of the *source transcript* (synthesis), not protein stock, is the supply variable; net sum over transcripts | supports |
| Yewdell 2003 [Y03] | presentation is dynamic and inefficient: < 0.1 % of peptides surviving degradation are presented | the source repertoire exceeds the detected ligandome ~3000-fold, as expected; expression is one variable among processing and binding | supports |
| Apcher 2011 [A11]; Apcher 2013 [A13] | PTC-carrying (NMD) mRNAs and pre-spliced RNAs in the nucleus produce class I peptides during the pioneer round, before decay; blocking cap-dependent translation suppresses full-length proteins but only partly affects antigenic peptides | low-stability transcripts are not excluded or down-weighted; bulk RNA-seq captures this supply, while the steady-state proteome under-represents it | supports |
| Nesvizhskii & Aebersold 2005 [NA05] | protein inference: peptides shared between proteins cannot be assigned unambiguously | reference proteomics loses or merges peptides of multicopy genes and paralogs, one reason to use RNA-seq | supports |
| Chion 2026 [Ch26] | proteomics still has limited dynamic range and coverage and many missing values for low-abundance species; DIA only partly solves this | proteomics is not used as the expression gold standard (Limitations) | consistent |
| Vendramin 2026 [V26] | active NMD restricts presentation of non-canonical, frameshift-derived peptides; SMG1 inhibition reveals them | non-canonical NMD-derived peptides are outside the MANE self model; steady-state TPM of NMD targets is low, so they are not over-counted | consistent, scoped |
| Pearson 2016 [P16] | 25,270 MAPs from 18 individuals; 41 % of expressed protein-coding genes give no MAPs; transcript and protein features predict MAP generation | expression is necessary but not sufficient; recall (97.7–98.1 %) is high while precision is low | supports |
| Kubiniok 2022 [K22] | MHC-I peptides come preferentially from actively transcribed genes and from evolutionarily hyperconserved proteins; HLA allotypes contribute unevenly | canonical (MANE) proteome; expression-weighted 9-mers; personal HLA typing | supports |
| Upadhya 2022 [U22] | the moderate mRNA–protein correlation partly reflects protein measurement error | proteomics is not used as a gold standard; a proteomics comparison is listed as a limitation | consistent, scoped |
| Madern 2025 [M25] | ribosomes on one mRNA cooperate through transient collisions, which reduces pausing; translation dynamics vary per mRNA | translation efficiency is not modelled (naive model) | consistent, scoped |

## Reference and transcript model

| paper | finding | pipeline step | verdict |
|---|---|---|---|
| Morales 2022 [Mo22] | MANE Select: one representative, well-expressed, conserved transcript per protein-coding gene, matched RefSeq/Ensembl | one MANE transcript per gene as pseudo-genome; one protein per gene | supports |
| Cheetham 2020 [C20]; Lin 2025 [L25] | pseudogenes are transcribed, some regulate their parents, and a few encode peptides | pseudogene-derived reads are treated as inflation of the parent (×funique only with pseudogene evidence); translated pseudogene peptides are not modelled | consistent, scoped. Some MHC pseudogenes keep the CDS but have lost the parent promoter/UTRs (HLA-DRB6, DPA2, DQB3), while HLA-H/J/K are near-complete copies that can be transcribed and derepressed in cancer ([pseudogene_problem](pseudogene_problem/README.md), [HLApseudo_MAFFT](../HLApseudo_MAFFT/README.md)). The correction touches 1.6–2.1 % of 9-mers and does not change MS AUC ([tests](tests/README.md)) |

## Mappability and multimapping

| paper | finding | pipeline step | verdict |
|---|---|---|---|
| Derrien 2012 [D12] | mappability varies between gene classes; discarding multimapped reads underestimates paralog expression | funique per transcript from simulated reads; fcM keeps multimappers | supports |
| Robert & Watson 2015 [RW15] | 12 quantification methods misestimate hundreds of disease-relevant genes in multicopy families; two-stage analysis proposed (unique, then multi-mapped groups) | fc (unique) + fcM (fractional); `shared = 1 − fc/fcM` decides the correction | supports |
| Deschamps-Francoeur 2020 [DF20] | `featureCounts -M --fraction` splits multireads evenly: every read counted once, but signal is diluted between active and inactive copies | fcM for totals; dilution handled by excluding alt/PAR copies and folding alt copies into the primary | supports |
| Dobin 2013 [Do13]; Li 2009 [Li09] | STAR alignment; SAM/BAM, MAPQ 255 = unique | alignment; funique = MAPQ-255 fraction | supports |

## HLA

| paper | finding | pipeline step | verdict |
|---|---|---|---|
| Aguiar 2019 [Ag19] | a personalised HLA index recovers more reads; reference-based estimates still correlate (r = 0.87, except HLA-DQA1) | personal MANE with the sample's HLA-A/B/C CDS; genome-wide counts unchanged, HLA −9 to +16 % | supports |
| Kang 2023 [Ka23] | personalised reference genomes are needed for accurate (single-cell) HLA expression | personal HLA records | supports |
| Boegel 2012 [B12]; Szolek 2014 [Sz14]; Orenbuch 2020 [O20] | HLA typing from RNA-seq; OptiType 97 % accuracy; arcasHLA 100 % at two fields for class I on the benchmark | OptiType + arcasHLA, reconciled. Lower agreement in our tests comes from 2 M-pair subsets and rare alleles, not tool failure | consistent, scoped |
| Barker 2023 [Ba23] (IPD-IMGT/HLA) | curated allele sequences | allele CDS from release 3.63.1 | supports |
| Bairoch 2018 [Ba18] (Cellosaurus) | cell-line knowledge resource that aggregates published HLA typings | tie-breaker only; conflicting records ranked by read support | supports |

## Data access

| paper | finding | pipeline step | verdict |
|---|---|---|---|
| Leinonen 2011 [Le11]; Maurya 2023 [Ma23] | the SRA archive; automated exploration of SRA datasets | `fetch` rule (prefetch / fasterq-dump) | supports |

**No paper contradicts a pipeline step.** Four findings mark scope limits of the naive self model:
non-canonical NMD peptides [V26], translated pseudogenes [C20, L25], translation efficiency [M25], and
mRNA–protein discordance [U22].

## References

- [Ag19] Aguiar VRC, César J, Delaneau O, Dermitzakis ET, Meyer D. Expression estimation and eQTL mapping for HLA genes with a personalized pipeline. *PLoS Genet.* 2019;15(4):e1008091. doi:10.1371/journal.pgen.1008091
- [A11] Apcher S, Daskalogianni C, Lejeune F, Manoury B, Imhoos G, Heslop L, Fåhraeus R. Major source of antigenic peptides for the MHC class I pathway is produced during the pioneer round of mRNA translation. *Proc Natl Acad Sci USA.* 2011;108(28):11572–11577. doi:10.1073/pnas.1104104108
- [A13] Apcher S, Millot G, Daskalogianni C, Scherl A, Manoury B, Fåhraeus R. Translation of pre-spliced RNAs in the nuclear compartment generates peptides for the MHC class I pathway. *Proc Natl Acad Sci USA.* 2013;110(44):17951–17956.
- [Ba18] Bairoch A. The Cellosaurus, a cell-line knowledge resource. *J Biomol Tech.* 2018;29(2):25–38.
- [Ba23] Barker DJ, Maccari G, Georgiou X, Cooper MA, Flicek P, Robinson J, Marsh SGE. The IPD-IMGT/HLA Database. *Nucleic Acids Res.* 2023;51(D1):D1053–D1060. doi:10.1093/nar/gkac1011
- [B12] Boegel S, Löwer M, Schäfer M, Bukur T, de Graaf J, Boisguérin V, Türeci Ö, Diken M, Castle JC, Sahin U. HLA typing from RNA-Seq sequence reads. *Genome Med.* 2012;4(12):102.
- [C20] Cheetham SW, Faulkner GJ, Dinger ME. Overcoming challenges and dogmas to understand the functions of pseudogenes. *Nat Rev Genet.* 2020;21:191–201. doi:10.1038/s41576-019-0196-1
- [Ch26] Chion M, Costanzo M, Gebelin M, et al. Proteomics at scale: bottlenecks and opportunities for early-career researchers in a fast developing field. *J Proteomics.* 2026;331:105679. doi:10.1016/j.jprot.2026.105679
- [D12] Derrien T, Estellé J, Marco Sola S, Knowles DG, Raineri E, Guigó R, Ribeca P. Fast computation and applications of genome mappability. *PLoS ONE.* 2012;7(1):e30377.
- [DF20] Deschamps-Francoeur G, Simoneau J, Scott MS. Handling multi-mapped reads in RNA-seq. *Comput Struct Biotechnol J.* 2020;18:1569–1576. doi:10.1016/j.csbj.2020.06.014
- [Do13] Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics.* 2013;29(1):15–21. doi:10.1093/bioinformatics/bts635
- [Ka23] Kang JB, Shen AZ, Gurajala S, et al. Mapping the dynamic genetic regulatory architecture of HLA genes at single-cell resolution. *Nat Genet.* 2023;55:2255–2268. doi:10.1038/s41588-023-01586-6
- [K22] Kubiniok P, Marcu A, Bichmann L, et al. Understanding the constitutive presentation of MHC class I immunopeptidomes in primary tissues. *iScience.* 2022;25(2):103768. doi:10.1016/j.isci.2022.103768
- [Le11] Leinonen R, Sugawara H, Shumway M. The Sequence Read Archive. *Nucleic Acids Res.* 2011;39:D19–D21. doi:10.1093/nar/gkq1019
- [Li09] Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. *Bioinformatics.* 2009;25(16):2078–2079. doi:10.1093/bioinformatics/btp352
- [L25] Lin YH, Yeh CT, Chen CY, Lin KH. Pseudogene: relevant or irrelevant? *Biomed J.* 2025;48:100790. doi:10.1016/j.bj.2024.100790
- [M25] Madern MF, Yang S, Witteveen O, Segeren HA, Bauer M, Tanenbaum ME. Long-term imaging of individual ribosomes reveals ribosome cooperativity in mRNA translation. *Cell.* 2025;188:1896–1911. doi:10.1016/j.cell.2025.01.016
- [Ma23] Maurya A, Szymanski M, Karlowski WM. ARA: a flexible pipeline for automated exploration of NCBI SRA datasets. *GigaScience.* 2023;12:giad067. doi:10.1093/gigascience/giad067
- [Mo22] Morales J, Pujar S, Loveland JE, et al. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research. *Nature.* 2022;604:310–315. doi:10.1038/s41586-022-04558-8
- [NA05] Nesvizhskii AI, Aebersold R. Interpretation of shotgun proteomic data: the protein inference problem. *Mol Cell Proteomics.* 2005;4(10):1419–1440.
- [O20] Orenbuch R, Filip I, Comito D, Shaman J, Pe'er I, Rabadan R. arcasHLA: high-resolution HLA typing from RNAseq. *Bioinformatics.* 2020;36(1):33–40. doi:10.1093/bioinformatics/btz474
- [P16] Pearson H, Daouda T, Granados DP, Durette C, Bonneil E, Courcelles M, Rodenbrock A, Laverdure JP, Côté C, Mader S, Lemieux S, Thibault P, Perreault C. MHC class I–associated peptides derive from selective regions of the human genome. *J Clin Invest.* 2016;126(12):4690–4701. doi:10.1172/JCI88590
- [RW15] Robert C, Watson M. Errors in RNA-Seq quantification affect genes of relevance to human disease. *Genome Biol.* 2015;16:177. doi:10.1186/s13059-015-0734-x
- [S00] Schubert U, Antón LC, Gibbs J, Norbury CC, Yewdell JW, Bennink JR. Rapid degradation of a large fraction of newly synthesized proteins by proteasomes. *Nature.* 2000;404(6779):770–774. doi:10.1038/35008096
- [Sz14] Szolek A, Schubert B, Mohr C, Sturm M, Feldhahn M, Kohlbacher O. OptiType: precision HLA typing from next-generation sequencing data. *Bioinformatics.* 2014;30(23):3310–3316. doi:10.1093/bioinformatics/btu548
- [U22] Upadhya SR, Ryan CJ. Experimental reproducibility limits the correlation between mRNA and protein abundances in tumor proteomic profiles. *Cell Rep Methods.* 2022;2:100288. doi:10.1016/j.crmeth.2022.100288
- [V26] Vendramin R, Fu H, Fernandez Patel S, Zhao Y, et al. Nonsense-mediated mRNA decay inhibition reshapes the cancer immunopeptidome. *Immunity.* 2026;59:1398–1421. doi:10.1016/j.immuni.2026.02.005
- [Y96] Yewdell JW, Antón LC, Bennink JR. Defective ribosomal products (DRiPs): a major source of antigenic peptides for MHC class I molecules? *J Immunol.* 1996;157(5):1823–1826.
- [Y03] Yewdell JW, Reits E, Neefjes J. Making sense of mass destruction: quantitating MHC class I antigen presentation. *Nat Rev Immunol.* 2003;3:952–961.
