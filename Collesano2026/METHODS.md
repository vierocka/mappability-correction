# Methods

Each step: what is done, why, where (script / rule) and the supporting reference.
References in **bold** are from the paper's reference list; *italic* ones are additional (not in the list).

## 1. Reads

| step | what | why | where | ref |
|---|---|---|---|---|
| download | `prefetch` + `fasterq-dump --split-files` (sra-tools 3.0.5); R1 ≠ R2 check | the 2026 runs passed R1 twice → HLA typing failed in 2 of 3 lines | `server/rules/qc.smk: fetch` | **Leinonen 2011** |
| QC | FastQC 0.12.1 text output parsed: mean/mode length per mate, read-weighted mean Phred (≥ 30 expected) | read length chooses the funique table; quality gate | `parse_fastqc.py`; `rule fastqc` | *Andrews 2010* |
| trimming (optional) | Trimmomatic 0.40 PE, ILLUMINACLIP TruSeq3-PE 2:30:10, LEADING 3, TRAILING 3, SLIDINGWINDOW 4:15, MINLEN 36 | adapter/quality filter; tested: < 4 % count change (`tests/`); paper aligned untrimmed reads (`align_trimmed: false`) | `rule trimmomatic` | *Bolger 2014* |

## 2. HLA typing

| step | what | why | where | ref |
|---|---|---|---|---|
| OptiType | razers3 fishing (`-i 95 -m 1 -dr 0`) against OptiType's RNA reference; ILP; `--rna` | genome-free: HLA reads are not lost to pseudogenes, alt contigs or patches | `rule optitype_fish`, `rule optitype` | **Szolek 2014**, *Weese 2012* |
| arcasHLA | `genotype -g A,B,C` from FASTQ; default `min_count` (75); reference rebuilt from IPD-IMGT/HLA 3.63.1; kallisto 0.44 is its internal engine (not used for quantification) | second, independent typer | `rule arcashla` | **Orenbuch 2020**, **Robinson 2019**, *Bray 2016* |
| reconcile | candidates = both typers + Cellosaurus; fished reads aligned to each candidate CDS (BWA-MEM); genotype support = reads explained; site check at MAFFT columns where candidates differ (allele base < 15 % median or < 5 % at a covered site → rejected); near tie (< 5 reads / 1 %) → best-supported Cellosaurus record; plot | the typers disagree on rare alleles (B\*07:412, C\*05:185, A\*24:08 in tests); Cellosaurus holds conflicting records | `reconcile_hla.py`; `rule reconcile_hla` | **Bairoch 2018**, *Katoh 2013*, *Li 2013* |
| group | runs of one group: calls of the run with most HLA reads; concordance table | technical replicates must agree | `group_hla_calls.py` | — |
| alleles | 2-field call → lowest-numbered full-length IMGT allele (CDS length = MANE CDS, ATG…stop) | one deterministic representative; identical protein | `select_hla_alleles.py` | **Robinson 2019** |
| personal MANE | MANE 5′/3′ UTR kept, HLA-A/B/C CDS replaced by the sample's alleles; heterozygous → two allele-suffixed records | allele-specific k-mers, correct peptides; UTR reads counted like every other gene | `build_personal_MANE_utr.py`, `build_personal_MANE_protein.py` | **Aguiar 2019**, **Morales 2022**, **Kang 2023** |

## 3. Quantification

| step | what | why | where | ref |
|---|---|---|---|---|
| reference | personal MANE Select v1.5 transcriptome as pseudo-genome (one transcript per gene); STAR index `--genomeSAindexNbases 12 --genomeChrBinNbits 11` | no alt/patch/pseudogene competitors (APM_locus: TAP1, PSMB8/9, HLA-DP = 0 on genome + alt) | `rule star_index_mane` | **Morales 2022** |
| alignment | STAR 2.7.10b `--alignIntronMax 1 --alignIntronMin 2 --outFilterMultimapNmax 20 --outSAMattributes NH HI AS NM MD --runRNGseed 0` | transcriptome alignment without splicing | `rule star_featurecounts` | **Dobin 2013** |
| counting | featureCounts 2.1.1 `-F SAF -p -B -C --countReadPairs -s 0`: fc (unique) and fcM (`-M --fraction`) | fcM keeps multimappers as fractions; fc gives the `shared` test and allele ratios | `rule star_featurecounts`; `make_mane_saf.py` | **Liao 2014**, **Deschamps-Francoeur 2020**, **Li 2009** |

## 4. Mappability correction

| step | what | why | where | ref |
|---|---|---|---|---|
| funique | MANE exons ± 50 bp tiled at 1-bp step, pseudo-pairs (R2 = revcomp R1), STAR 2.7.11b vs Ensembl GRCh38 primary assembly (194 seqs) `--outFilterMultimapNmax 40 --outFilterMismatchNmax 1 --outSAMmultNmax 1`; funique = MAPQ-255 fraction of R1; tables L75 / L100, other lengths simulated | sequence uniqueness per transcript at the library's read length | `compute_uniqueness.py`, `build_correction_tables.py`; `funique.smk` | **Derrien 2012**, **Storvall 2013** |
| exclusions | 70 non-coding, 17 PAR, 7 alt copies folded into the primary, 20 alt-only population-variable genes (HLA-DRB3/4 kept for DR52/DR53), 3 without a unique locus | no peptides / not quantifiable per copy / not in primary assembly | `excluded_transcripts/`; `apply_uniqueness_and_average.py` | **Morales 2022**, **Robert & Watson 2015** |
| ×funique | only if f < 0.9 **and** pseudogene evidence **and** shared = 1 − fc/fcM < 0.10 (or fcM < 20); f = 0 → 0.001 | inflation by reads from paralogs/pseudogenes outside MANE | same | **Robert & Watson 2015**, **Deschamps-Francoeur 2020** |
| HLA | total = Σ allele fcM; split by unique-read ratio | heterozygous alleles share reads; unique-only loses 19–82 % | same (`hla_total_source`) | **Aguiar 2019** |
| TPM | counts / length (kb), scaled to 10⁶ over retained transcripts; replicate mean per condition | comparable abundance | same | *Wagner 2012* |

## 5. Peptides (local)

| step | what | why | where | ref |
|---|---|---|---|---|
| 9-mers | all 9-mers of the personal MANE protein per transcript (HLA allele IDs kept) | candidate class I ligands | `generate_maneSelect_9mer_table_1gene_1canonicalProtein.py`, `ninemer_table_to_rows.py` | **Morales 2022** |
| net sum | Σ transcript TPM_corrected × occurrences; transcripts with TPM_corrected < 1 dropped; peptides with U/X dropped | peptide supply proxy | `rebuild_peptide_tables.py` | — |
| reports | loss per filter step (% counts, % TPM), composition shift; inflation tail; recall of MS ligands (sanity check only — saturated near its ceiling) and ROC-AUC of net-sum ranks | interpretability | `report_filter_steps.py`, `inflation_tail.py`, `impact_filters.py`, `compare_ms_bassani.py` | *Bassani-Sternberg 2015* |

## References

From the paper's list:
- **Leinonen R**, Sugawara H, Shumway M. The Sequence Read Archive. Nucleic Acids Res 2011;39:D19–21.
- **Deschamps-Francoeur G**, Simoneau J, Scott MS. Handling multi-mapped reads in RNA-seq. Comput Struct Biotechnol J 2020;18:1569–76.
- **Robert C**, Watson M. Errors in RNA-Seq quantification affect genes of relevance to human disease. Genome Biol 2015;16:177.
- **Kang JB** et al. Mapping the dynamic genetic regulatory architecture of HLA genes at single-cell resolution. Nat Genet 2023;55:2255–68.
- **Derrien T** et al. Fast computation and applications of genome mappability. PLoS ONE 2012;7:e30377.
- **Storvall H**, Ramsköld D, Sandberg R. Efficient and comprehensive representation of uniqueness for next-generation sequencing by minimum unique length analyses. PLoS ONE 2013;8:e53822.
- **Aguiar VRC** et al. Expression estimation and eQTL mapping for HLA genes with a personalized pipeline. PLoS Genet 2019;15:e1008091.
- **Boegel S** et al. HLA typing from RNA-Seq sequence reads. Genome Med 2012;4(12):102 (listed as 2013 in the paper's list; the volume is 2012).
- **Szolek A** et al. OptiType: precision HLA typing from next-generation sequencing data. Bioinformatics 2014;30:3310–6.
- **Orenbuch R** et al. arcasHLA: high-resolution HLA typing from RNAseq. Bioinformatics 2020;36:33–40.
- **Bairoch A**. The Cellosaurus, a cell-line knowledge resource. J Biomol Tech 2018;29:25–38.
- **Robinson J** et al. IPD-IMGT/HLA Database. Nucleic Acids Res 2019.
- *Barker DJ* et al. The IPD-IMGT/HLA Database. Nucleic Acids Res 2023;51:D1053–60 (current database paper).
- **Morales J** et al. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research. Nature 2022;604:310–5.
- **Dobin A** et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 2013;29:15–21.
- **Liao Y**, Smyth GK, Shi W. featureCounts. Bioinformatics 2014;30:923–30.
- **Li H** et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 2009;25:2078–9.

Additional (not in the paper's list; add if the step is kept). The biological rationale (DRiPs, pioneer round, NMD,
immunopeptidome determinants, pseudogenes) is cited in [LITERATURE.md](LITERATURE.md).
- *Andrews S*. FastQC. Babraham Bioinformatics, 2010.
- *Bolger AM*, Lohse M, Usadel B. Trimmomatic. Bioinformatics 2014;30:2114–20.
- *Weese D*, Holtgrewe M, Reinert K. RazerS 3. Bioinformatics 2012;28:2592–9.
- *Bray NL* et al. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol 2016;34:525–7 (kallisto, inside arcasHLA).
- *Katoh K*, Standley DM. MAFFT version 7. Mol Biol Evol 2013;30:772–80.
- *Li H*. Aligning sequence reads … with BWA-MEM. arXiv:1303.3997, 2013.
- *Wagner GP*, Kin K, Lynch VJ. Measurement of mRNA abundance using RNA-seq data. Theory Biosci 2012;131:281–5.
- *Bassani-Sternberg M* et al. Mass spectrometry of HLA class I peptidomes … Mol Cell Proteomics 2015;14:658–73.
- *Mölder F* et al. Sustainable data analysis with Snakemake. F1000Research 2021;10:33.
- *Ewels PA* et al. The nf-core framework. Nat Biotechnol 2020;38:276–8 (primary-assembly recommendation).
