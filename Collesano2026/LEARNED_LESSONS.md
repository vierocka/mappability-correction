# Lessons: why the final pipeline looks like this

This file condenses the earlier tests behind the final pipeline. The literature check is in [LITERATURE.md](LITERATURE.md). Every point links to the place where it was investigated.
Tests from project folders outside this repository (`human/RNAseq`, `human/References`) were copied, with paths
sanitised, to [earlier_tests/](earlier_tests/README.md).

## 1. The self model: deliberately naive

**Self = filtered MANE Select (v1.5): one canonical protein per gene, with personal HLA-A/B/C.**

| assumption | why | investigated in |
|---|---|---|
| one protein per gene (MANE Select) | MANE: 19,289 genes, 10.94 M unique 9-mers. All RefSeq isoforms: 19,866 genes, 12.09 M (+10 %); median 421 vs 458 9-mers per gene. Extra isoforms add few new 9-mers but split the expression. | [peptides_per_gene_MANE_vs_GRCh38p14.py](earlier_tests/References_peptides/peptides_per_gene_MANE_vs_GRCh38p14.py), [notes](earlier_tests/References_peptides/UniquePeptides_per_gene_MAN_vs_GRCh38p14_notes.txt), [Suppl. Fig. 1B](results/figures/isoforms/SupplFig1_isoforms_ABCD.png) |
| expression = net sum over transcripts | a 9-mer is supplied by every transcript that encodes it (Σ TPM × occurrences) | [rebuild_peptide_tables.py](scripts/rebuild_peptide_tables.py), [METHODS §5](METHODS.md#5-peptides-local) |
| keep low-stability and decay-directed transcripts | defective ribosomal products and transcripts translated before decay (pioneer round, NMD substrates, pre-spliced RNA) are a major source of class I peptides, and presentation is gene-selective; mRNA abundance is not protein abundance [1–5]. Active NMD still limits presentation of *non-canonical* (frameshift) peptides (Vendramin 2026); these are outside the MANE self model | [References](#references), [LITERATURE](LITERATURE.md) |
| no protein-level filter | a proteomics comparison is still to do | [README: limitations](../README.md#limitations) |
| not too inflated | pseudogene/paralog reads corrected (×funique), non-coding and multicopy transcripts removed, transcript TPM ≥ 1. Pseudogenes can be transcribed and a few translated (Cheetham 2020; Lin 2025). Some keep the CDS but have lost the parent promoter/UTRs (HLA-DRB6, DPA2, DQB3), so they add few reads. Class I HLA-H/J/K are near-complete copies that keep promoter and 3′UTR ([pseudogene_problem](pseudogene_problem/README.md), [HLApseudo_MAFFT](../HLApseudo_MAFFT/README.md)). Their own peptides are not modelled | [excluded_transcripts](excluded_transcripts/README.md), [apply_uniqueness_and_average.py](scripts/apply_uniqueness_and_average.py), [tests/filter_impact](tests/README.md), [tests/inflation_tail](tests/README.md) |
| still self | personal HLA alleles. Somatic variants are not included: the HCC1937 nonsense/missense gene classes were listed but not used | [select_hla_alleles.py](scripts/select_hla_alleles.py), [gene_classes_notes.txt](earlier_tests/RNAseq_isoforms/gene_classes_notes.txt) |

**Not modelled:** gene-specific degradation (RNA integrity, 3′ bias), §4.

<a id="isoforms"></a>
## 2. Isoforms: kallisto vs STAR + featureCounts

| test (HCC1937 SRR925709 unless stated) | result | investigated in |
|---|---|---|
| kallisto 0.50.1 on RefSeq GRCh38.p14 CDS (146,337 targets, k = 31) | 60.3 % of 24.6 M pairs pseudoaligned | [hpc_kallisto.sh](earlier_tests/RNAseq_isoforms/hpc_kallisto.sh), [mapping_overview](earlier_tests/RNAseq_isoforms/HCC1937_mapping_overview.txt) |
| STAR 2.7.11b on the full GRCh38.p14 genome (with alt/patches) | 80.0 % unique, 15.6 % multi | [STAR Log.final.out](earlier_tests/RNAseq_isoforms/HCC1937_STAR_Log.final.out) |
| featureCounts, unique only vs `-M` | 41.5 % vs 83 % of reads assigned | [featureCounts summary](earlier_tests/RNAseq_isoforms/HCC1937_featureCounts.summary), [mapping_overview](earlier_tests/RNAseq_isoforms/HCC1937_mapping_overview.txt) |
| patch/alt duplicates in the CDS index | 14,069 targets share an accession (max taken); `seqkit rmdup -s` removed 53,249 identical CDS | [dedup_kallisto_accessions.py](earlier_tests/RNAseq_isoforms/dedup_kallisto_accessions.py), [deduplication note](earlier_tests/References_peptides/GRCh38p14_deduplication_note.txt) |
| isoforms per expressed gene | ≥ 2: 59 %, ≥ 5: 15 %, ≥ 10: 3 % | [JY_HCC1937_9merRepertoires.R](earlier_tests/RNAseq_isoforms/JY_HCC1937_9merRepertoires.R), [table](results/figures/isoforms/HCC1937_isoforms_per_gene_fraction_ge_k.tsv) |
| top isoform share (kallisto, ≥ 10 counts) | JY / HCC1937: only 34–35 % of genes have max/sum ≥ 0.9; 36 % < 0.5 | [gene_isoform_sum_max.py](earlier_tests/RNAseq_isoforms/gene_isoform_sum_max.py), [max_vs_sum.R](earlier_tests/RNAseq_isoforms/max_vs_sum.R), [plot_makie_sum_vs_max_kallisto.jl](earlier_tests/RNAseq_isoforms/plot_makie_sum_vs_max_kallisto.jl), [JY](results/figures/isoforms/JY_sum_vs_max_makie.png) · [HCC1937](results/figures/isoforms/HCC1937_sum_vs_max_makie.png) · [HCT116](results/figures/isoforms/HCT116_sum_vs_max_makie.png) |
| **ties** | 31 % of multi-isoform genes have the top two isoforms within 1 %; spike at exactly 0.50 in two-isoform genes (≈ 600 genes) | [isoform_ties.py](earlier_tests/RNAseq_isoforms/isoform_ties.py), [JY_HCC1937_9merRepertoires.R](earlier_tests/RNAseq_isoforms/JY_HCC1937_9merRepertoires.R), [plot](results/figures/isoforms/HCC1937_two_isoform_major_fraction.png) |
| kallisto max isoform vs featureCounts gene | Spearman 0.81; median ratio 0.26; 735 genes 0 in fc but ≥ 10 in kallisto, 52 the reverse | [max_vs_sum.R](earlier_tests/RNAseq_isoforms/max_vs_sum.R), [merged table](earlier_tests/RNAseq_isoforms/fc_vs_kallisto_max_merged.tsv), [isoform_ties.py](earlier_tests/RNAseq_isoforms/isoform_ties.py), [scatter](results/figures/isoforms/fc_vs_kallisto_max_scatter.png) |
| unique 9-mers, canonical (fc) vs leading isoform (kallisto) | 9,200 genes same, 7,201 more, 3,675 fewer | [canonicalFCset_leadingIsoformKallisto_plots.R](earlier_tests/RNAseq_isoforms/canonicalFCset_leadingIsoformKallisto_plots.R), [table](earlier_tests/RNAseq_isoforms/HCC1937_unique9mers_per_gene_canonical_vs_leading_table.tsv) |
| isoforms on simulated reads (dedup-CDS index) | 19,423 transcripts: 47 % > 10 % too low, 21 % zero, 6 % > 10 % too high; TAPBP (14 isoforms) −44 % | [01_isoform_redundancy](../mappability/smallPhenomena/01_isoform_redundancy/), [APM_locus §3](APM_locus/README.md) |

**Reading.** The apparent multi-isoform expression is mostly EM ambiguity: identical or near-identical CDS cannot be
told apart, so kallisto splits them (the exact 0.50 spike). Few genes need more than one isoform. kallisto spreads a
gene into fractions; STAR + featureCounts with unique counting gives zeros wherever copies exist.
**Choice:** one MANE transcript per gene, STAR with `fcM` fractions, and a correction only where a pseudogene explains
low uniqueness.

## 3. UTRs matter in bulk RNA-seq

| test (simulated reads) | result | investigated in |
|---|---|---|
| CDS reads → CDS index | 5–6 % unique | [kallisto_results](../mappability/kallisto_results/README.md) |
| CDS reads → MANE (with UTRs) | 95–96 % unique | [kallisto_results](../mappability/kallisto_results/README.md) |
| unique k-mers, CDS → MANE RNA | 3 % → 94 % | [Why UTRs fail…](../mappability/smallPhenomena/01_isoform_redundancy/Why_UTRs_fail_to_rescue_nearly2000annotGenes.txt) |
| genes not rescued by UTRs | 1,922: large families with slow UTR drift (OR, ZNF, KRTAP, PCDH, KIR, ≈ 1,600); en-bloc duplications (PRR20, DEFB, AMY, ≈ 200); histones without a polyA 3′UTR (≈ 40); Pol III ncRNA (RNU6, ≈ 17) | [Why UTRs fail…](../mappability/smallPhenomena/01_isoform_redundancy/Why_UTRs_fail_to_rescue_nearly2000annotGenes.txt), [RNU_snRNA](../RNU_snRNA/README.md) |
| side effect | 71 % of transcripts are "MANE-inflated": UTR k-mers rescue reads that are ambiguous in the CDS. 5 % are CDS-inflated, because their UTRs are shared (HBA1/HBA2) | [kallisto_results](../mappability/kallisto_results/README.md) |
| 3′-end protocols (scRNA) | uniqueness differs when only the 3′ end is sequenced | [sc3UTRs](../mappability/sc3UTRs/README.md) |

**Choice:** align to full MANE transcripts (UTRs kept, also in the personal HLA records); translate the CDS only.

## 4. Degradation: out of scope, but it matters

Partly degraded RNA gives 3′-biased coverage. A truncated transcript can look like a shorter isoform, and EM then
moves reads between isoforms. This adds to the tie problem in §2: an exact 0.50 split points to identical sequence,
while skewed splits in degraded libraries can come from 3′ bias. With one MANE transcript per gene the effect
shrinks to gene level (lower TPM for long genes in degraded samples).
*Investigated in a separate degradation-QC project, not in this repository.*

## 5. HLA and the APM locus

| point | investigated in |
|---|---|
| Genome-based HLA typing lost HLA reads to pseudogenes, alt haplotypes and patches (and R1 was passed twice); typing is now genome-free (OptiType, arcasHLA) | [METHODS §2](METHODS.md#2-hla-typing), [original_shell](scripts/original_shell/) |
| typers disagree on rare alleles (JY B\*07:412 vs B\*07:02); Cellosaurus disagrees with itself (HCC1937 B\*40:16/81:03); reconcile uses read support and base checks, with Cellosaurus only as tie-breaker | [reconcile_hla.py](scripts/reconcile_hla.py), [tests/hla_typing](tests/hla_typing/), [tests README](tests/README.md) |
| HLA pseudogenes share most 100-bp windows with HLA-A (63 %) and DRB1 (66 %) | [HLApseudo_MAFFT](../HLApseudo_MAFFT/README.md), [HLA_MAFFT](../HLA_MAFFT/README.md), [pseudogene decoys](pseudogene_problem/decoys.md) |
| APM locus: alt haplotypes give zeros, antisense/nested genes overlap, isoforms split | [APM_locus](APM_locus/README.md), [Figure_APMlocus](../Figure_APMlocus/README.md) |

## 6. Other observations

| point | investigated in |
|---|---|
| HCC1937 is breast cancer, HCT116 colon cancer (an old cluster script had the two labels swapped in a comment) | [cancer_cell_lines.csv](earlier_tests/RNAseq_isoforms/cancer_cell_lines.csv), [hpc_kallisto.sh](earlier_tests/RNAseq_isoforms/hpc_kallisto.sh) |
| the earlier STAR run used reads trimmed to 2 × 70 bp (index sjdbOverhang 70); the final counts use untrimmed reads | [STAR Log.final.out](earlier_tests/RNAseq_isoforms/HCC1937_STAR_Log.final.out), [tests/trim_personal](tests/README.md) |
| the peptide sets of the cell lines differ mainly in expression, not sequence; only HLA is personalised | [build_personal_MANE_utr.py](scripts/build_personal_MANE_utr.py) |

## References

1. Yewdell JW, Antón LC, Bennink JR. Defective ribosomal products (DRiPs): a major source of antigenic peptides for
   MHC class I molecules? *J Immunol.* 1996;157(5):1823–1826.
2. Schubert U, Antón LC, Gibbs J, Norbury CC, Yewdell JW, Bennink JR. Rapid degradation of a large fraction of newly
   synthesized proteins by proteasomes. *Nature.* 2000;404(6779):770–774. doi:10.1038/35008096
3. Apcher S, Daskalogianni C, Lejeune F, Manoury B, Imhoos G, Heslop L, Fåhraeus R. Major source of antigenic
   peptides for the MHC class I pathway is produced during the pioneer round of mRNA translation.
   *Proc Natl Acad Sci USA.* 2011;108(28):11572–11577. doi:10.1073/pnas.1104104108
4. Apcher S, et al. Translation of pre-spliced RNAs in the nuclear compartment generates peptides for the MHC class I
   pathway. *Proc Natl Acad Sci USA.* 2013;110(44):17951–17956. PMCID: PMC3816431
5. Pearson H, et al. MHC class I–associated peptides derive from selective regions of the human genome.
   *J Clin Invest.* 2016;126(12):4690–4701. doi:10.1172/JCI88590
