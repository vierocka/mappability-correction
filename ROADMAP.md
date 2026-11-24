# ROADMAP: mappability-correction (file guide)

One entry per file. Repeated families use a `*` / `{a..o}` pattern.

## Top level
| file | what |
|---|---|
| `README.md` | Project intro: in-silico mappability correction for RNA-seq at the HLA/APM locus; repo layout. |
| `README (… conflicted copy 2026-06-18).md` | Dropbox conflict copy of the README (older, longer). Merge or delete. |
| `ROADMAP.md` | This file guide. |
| `comparison_scheme.md` | What is compared and why; escalation path when short reads are insufficient. |
| `deduplication_note.txt` | `seqkit rmdup -s` command used to make the deduplicated RefSeq CDS FASTA. |
| `shortRNAseq_consortia_projects.txt` | Table of large short-read RNA-seq consortia (TCGA …) that the correction applies to. |
| `LICENSE` | Repository licence. |
| `.gitignore`, `.gitignore (… conflicted copy)` | Git ignore rules; the second is a Dropbox conflict copy. |
| `versions/versions.txt` | HPC module versions (STAR 2.7.11b, SAMtools 1.19.2, Python 3.11.5). |

## `Collesano2026/`: reproducible pipeline of the paper
| file | what |
|---|---|
| `README.md`, `LEARNED_LESSONS.md`, `METHODS.md`, `ROADMAP.md`, `PHENOMENA.md` | Overview + validation; every method with citation; decisions and known discrepancies; phenomena table + figure. |
| `config/server.yaml`, `config/samples.tsv`, `config/local.yaml` | Pipeline settings; one row per run (SRA or FASTQ, group, Cellosaurus ID); local inputs. |
| `server/Snakefile`, `server/rules/{qc,align,correct,report}.smk`, `server/helpers/download_references.sh` | Server pipeline: download/QC/trim, STAR + featureCounts, correction, summaries; reference download. |
| `local/Snakefile` | Local pipeline: 9-mers, net sums (TPM ≥ 1), filter report, inflation tail, MS comparison. |
| `scripts/group_hla_calls.py`, `scripts/summarize_server.py`, `scripts/make_mane_saf.py`, `scripts/plot_phenomena.py` | Group HLA calls; server summary tables/plot; MANE SAF; phenomena figure. |
| `results/figures/phenomena_overview.*`, `results/tables/phenomena_*.tsv` | Phenomena overview figure and its tables. |
| `results/figures/isoforms/*` | Isoform evidence copied from the earlier RNAseq tests (Suppl. Fig. 1, sum vs max, ties, fc vs kallisto). |
| `earlier_tests/*` | Copied scripts, logs and tables of the earlier kallisto/STAR/featureCounts isoform tests and MANE vs RefSeq peptide sets (sanitised). |
| `docs/subway/*` | slitherFlow subway-map specs (yaml) and renders (svg, png) of both pipelines. |
| `excluded_transcripts/README.md` | The 117/115 excluded MANE transcripts: groups, reasons, references. |
| `excluded_transcripts/1_non_coding.tsv` | 70 non-coding MANE transcripts (no protein → no peptides). |
| `excluded_transcripts/2_pseudoautosomal_PAR.tsv` | 17 PAR1/PAR2 genes (identical X/Y copies, funique = 0). |
| `excluded_transcripts/3_alt_copies_folded_into_primary.tsv` | 7 alt-haplotype copies; their counts are added to the primary transcript. |
| `excluded_transcripts/4_population_variable_alt_haplotype.tsv` | 20 alt-only, population-variable genes (KIR …); HLA-DRB3/4 kept in JY. |
| `excluded_transcripts/5_no_unique_primary_locus.tsv` | 3 fix-patch transcripts without a unique primary locus. |
| `excluded_transcripts/6_reincluded_non_primary_f1.tsv` | 34 non-primary transcripts kept with f = 1. |
| `resources/non_primary_MANE_transcripts.tsv` | Curated list of non-primary MANE transcripts (`excluded`, `notes`); verbatim copy. |
| `resources/non_coding_MANE_transcripts_no_protein.tsv` | MANE transcripts without protein; verbatim copy. |
| `scripts/PROVENANCE.tsv` | Origin path + md5 of every vendored script. |
| `scripts/apply_uniqueness_and_average.py` | Original correction (exclusions, fold, ×funique, TPM, replicate mean); env hooks added. |
| `scripts/make_hla_tap_tables.py` | Original HLA/TAP summary tables; config hook added. |
| `scripts/rebuild_peptide_tables.py` | Original 9-mer net-sum tables; hook adds transcript TPM ≥ 1 cut-off. |
| `scripts/build_personal_MANE_utr.py` | Original: personal MANE (MANE UTRs kept, HLA CDS replaced by IMGT alleles). |
| `scripts/build_personal_MANE_protein.py` | Original: personal MANE protein FASTA. |
| `scripts/generate_maneSelect_9mer_table_1gene_1canonicalProtein.py` | Original 9-mer generator; regex fixed to keep HLA allele suffixes. |
| `scripts/build_correction_tables.py`, `scripts/compute_uniqueness.py` | Original funique table builders (simulation → MAPQ-255 fraction). |
| `scripts/original_shell/*.sh` | Original HPC shell scripts (Laura personal MANE, arcasHLA, OptiType, mappability master). |
| `scripts/ninemer_table_to_rows.py` | New: expands the 9-mer table into one row per occurrence. |
| `scripts/parse_fastqc.py` | New: FastQC text → mean/mode read length, mean Phred (≥ 30 check), STAR length. |
| `scripts/report_filter_steps.py` | New: % counts/TPM lost per filter step; composition shift plots. |
| `scripts/compare_ms_bassani.py` | New: source-repertoire 9-mers vs Bassani-Sternberg 2015 MS (recall, AUC, fate of each ligand). |
| `scripts/impact_filters.py` | New: raw / exclusions / ×funique / TPM ≥ 1 variants scored against MS. |
| `scripts/inflation_tail.py` | New: net-sum histograms before/after funique; is the extreme tail corrected? |
| `envs/*.yaml` | Pinned conda envs: qc, align (STAR 2.7.10b), genome (STAR 2.7.11b), arcashla, optitype, py. |
| `APM_locus/README.md` | Why chr6 HLA/APM needs special handling: pseudogenes, alt contigs, overlaps, isoforms; STAR zero vs kallisto fractions; pipeline choices; open issues. |
| `APM_locus/star_funique_by_reference_L100.tsv`, `kallisto_truth_vs_dedupCDS.tsv` | APM-gene funique per reference; kallisto TPM, MANE vs dedup-CDS index. |
| `scripts/reconcile_hla.py` | OptiType vs arcasHLA: read support, site check (MAFFT, base quality), Cellosaurus; plot. |
| `scripts/select_hla_alleles.py` | 2-field calls → IMGT CDS FASTA (lowest full-length allele) for the personal MANE. |
| `scripts/choose_read_length.py` | Mean read length → funique L75 / L100 / simulate new length. |
| `resources/funique/*` | Precomputed funique tables L75, L100 (HPC 2026-08-12), plus cache for new lengths. |
| `server/rules/{funique,hla}.smk` | Snakemake rules: funique per read length; genome-free HLA typing → personal MANE. |
| `tests/*` | Impact tests: filtering, inflation tail, trimming, personal MANE, HLA typing, FastQC (2 M-pair subsets). |
| `server/`, `local/` | Snakemake pipelines (server: FASTQ → HLA + corrected TPM; local: TPM → 9-mer net sums). In progress. |

## `mappability/`: funique simulation study
| file | what |
|---|---|
| `README2.md` | Mappability study description (designs v1–v8, outputs). |
| `roadmap.csv` | One row per script version: read length, reference, read source, flanks, SE/PE, output TSV. |
| `scripts_overview.csv` | Shorter version table of v1–v5 designs; flags output-dir conflicts. |
| `check_the_path_to_your_result_per_script.tsv` | Lookup: script version → results path. |
| `reference_differences.docx` | Notes on differences between the reference builds (RefSeq, Ensembl, MANE, GENCODE). |
| `build_dedup_cds_nf.py` | Builds the deduplicated CDS FASTA on the Ensembl primary assembly ("nf" = nf-core-style reference). |
| `build_id_conversion_table.py` | RefSeq ↔ Ensembl ↔ MANE ID conversion table. |
| `check_exon_positions_nf.py`, `check_id_conversion_nf.py` | QC of exon coordinates and ID mapping for the nf references. |
| `compare_cds_counts.py` | Compares CDS counts between reference builds. |
| `generate_scripts.py` | Generates the `Map_scripts/mappability_correction_v*.sh` family. |
| `summarise_results.py`, `summarise_results_cds.py` | Summarise transcript-/CDS-level uniqueness factors into `results/results_summary.tsv`, `cds_summary.tsv`. |
| `summarise_kallisto.py`, `summarise_kallisto_genes.py` | Summarise kallisto cross-validation runs. |
| `extract_dedup.log` | Log of the dedup-CDS extraction. |
| `kallisto_results.tar.gz` | Archive of `kallisto_results/`. |

### `mappability/Map_scripts/`
| file | what |
|---|---|
| `mappability_correction_v{1..8}_{a..o}.sh` | One simulation design each (see `roadmap.csv`): reads → STAR → funique TSV. |
| `mappability_correction_v*_nf.sh` | Same designs on the Ensembl GRCh38.116 primary assembly (used for the paper: v1_a_nf). |
| `mappability_v7_fix_uf.sh` | Fix-up of v7 uniqueness factors. |
| `mappability_v8_HPC.sh` | v8 ("proper PE") HPC submission script. |
| `run_all_mappability.sh`, `run_all_mappability_nf.sh` | Submit all versions (original / nf). |
| `kallisto_mappability_crossval.sh` | kallisto cross-validation of the simulated reads. |
| `generate_scripts_nf.py` | Generator for the `_nf` scripts. |
| `overview_v8_notes.txt`, `README3.md` | Notes on v8 and on the script family. |
| `__pycache__/` | Python bytecode cache (ignore). |

### `mappability/results/`
| file | what |
|---|---|
| `transcript_uniqueness_factors_<design>_L<len>bp[_nf\|_proper].tsv` | funique per MANE transcript for each design/read length (L75 = 75 bp, unlabeled = 100 bp). |
| `cds_uniqueness_factors_{genome,selfmap}_{SE,PE}[_L*].tsv` | funique per deduplicated CDS (v7). |
| `results_summary.tsv`, `cds_summary.tsv`, `cds_counts_by_gene.tsv` | Cross-design summaries. |
| `id_conversion_table.tsv`, `id_conversion_check_nf.tsv`, `exon_position_check_nf.tsv`, `chr_name_map.tsv` | ID and coordinate checks. |
| `mappability_out.logs`, `README.md` | Run logs; results description. |

### `mappability/kallisto_results/`
| file | what |
|---|---|
| `mappability_<design>_{SE,PE}[_L*]/kallisto_vs_{MANE,dedup_cds}/` | kallisto quantification of simulated reads per design. |
| `kallisto_run_summary.tsv`, `kallisto_run_overview.tsv`, `kallisto_genes_summary.tsv`, `kallisto_inflation.tsv` | Summaries: per run, per gene, inflation vs truth. |
| `chr_name_map.tsv`, `README.md` | Chromosome naming; description. |

### `mappability/sc3UTRs/`: 3′-end (single-cell) uniqueness
| file | what |
|---|---|
| `compute_3prime_uf.py`, `compute_3end_similarity.py` | funique restricted to 3′ ends; 3′-end similarity between genes. |
| `annotate_gene_overlaps.py`, `intersect_markers.py` | Gene overlap annotation; intersect with cell-type markers. |
| `download_markers.sh`, `CellMarker2_Human_cell_markers.xlsx`, `PanglaoDB_markers_27Mar2020.tsv` | Marker databases and download script. |
| `run_3prime_HPC.sh`, `README.md` | HPC job; description. |
| `results/sc3prime_v*_*.tsv` | 3′ funique per design; `celltype_impact_summary.tsv`, `marker_*.tsv` = marker impact. |

### `mappability/smallPhenomena/`: case studies
| file | what |
|---|---|
| `README.md`, `case_study_lib.py`, `lib_check.py` | Overview; shared helper library and its check. |
| `01_isoform_redundancy/*` | Isoform homology/k-mer/mappability tables and gene/cluster recommendations. |
| `02_overlapping_annotation/`, `04_multicopy_coding/`, `05_local_duplications/`, `06_unique_genes/`, `07_sparse_celltypes/` | `analyze.py` only (planned case studies). |
| `03_repetitive_pseudogene/*` | Pseudogene genes/summary, UTR benefit; reference PDF (gkn758). |
| `08A_immune_evasion/*`, `08B_HLA_MHC/*` | Gene classification, STAR/kallisto funique across runs, critique notes, article drafts. |

## `Ref/`: references (large, not in git)
| file | what |
|---|---|
| `Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa(.fai)` | Ensembl primary assembly (194 seqs); funique alignment target. |
| `Homo_sapiens.GRCh38.116.gtf` | Ensembl 116 annotation. |
| `Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds{,_mane}.fna` | Deduplicated CDS on the primary assembly (all / MANE). |
| `GCF_000001405.40_GRCh38.p14_genomic.{fna,fai,gff}` | RefSeq GRCh38.p14 full genome + alt/patches. |
| `GCF_000001405.40_GRCh38.p14_primary_only.fna(.fai)` | RefSeq primary chromosomes only. |
| `GCF_000001405.40_GRCh38.p14_cds_from_genomic.dedup.fna` | RefSeq CDS, deduplicated (see `deduplication_note.txt`). |
| `MANE.GRCh38.v1.5.ensembl_{genomic.gff,rna.fna}` | MANE Select v1.5 annotation and transcript FASTA. |
| `gencode.v49.basic.annotation.gtf` | GENCODE v49 basic annotation. |
| `recommendations_nf.txt` | Why the Ensembl primary assembly is used (nf-core recommendation). |

## Figure and alignment folders
| file | what |
|---|---|
| `Figure_APMlocus/plot_APM_locus.py` → `APM_locus_chr6.{pdf,png}` | APM locus map on chr6; `APM_locus_GFF_snapshot.tsv` = input; `README.md`. |
| `HLA_MAFFT/extract_HLA_sequences.py`, `run_mafft.sh`, `plot_alignments.py` | HLA sequence extraction (MANE / primary chr / dedup CDS), MAFFT, plots. |
| `HLA_MAFFT/HLA_{MANE,primaryChr,dedup_cds}.{fna,mafft.fna,mafft.log}` | Inputs, alignments and logs per source. |
| `HLA_MAFFT/HLA_alignments{,_A,_B}.{pdf,png}` | Alignment figures; `plot_notes.txt`, `TransparentScience_HLA_Li.txt`, `README.md`. |
| `HLApseudo_MAFFT/run_hla_mafft.py`, `shared_fraction.py`, `plot_shared_fraction.py` | HLA genes vs pseudogenes per clade: MAFFT, shared-position fraction, plot. |
| `HLApseudo_MAFFT/class{I_A,I_BC,II_DP,II_DQ,II_DRB}*.{input.fa,aln.fa,aln.clustal,aln.log}` | Per-clade inputs and alignments. |
| `HLApseudo_MAFFT/shared_fraction.tsv`, `shared_positions.tsv`, `transcript_table.tsv`, `hla_shared_fraction.{pdf,png}` | Results and figure; notes in `notes.txt`, `figure_notes.txt`, `explore_commands.sh`, `README.md`. |
| `RNU_snRNA/extract_RNU_sequences.py`, `run_mafft.sh`, `plot_RNU_homology.py`, `RNU_UF_summary.py` | RNU snRNA families: extraction, MAFFT, homology plot, funique summary. |
| `RNU_snRNA/RNU{1,2,4,5,6}_family.{fna,mafft.fna,mafft.log}`, `RNU_MANE{,_unique}.fna` | Family sequences/alignments; MANE RNU set. |
| `RNU_snRNA/RNU_homology.{pdf,png}`, `RNU_UF_summary.tsv`, `README.md` | Figure, funique table, description. |
| `Figure_MHC_homology/combine_panels_ABC.py` → `Figure_ABC{,_vchrom}.{pdf,png}` | Supplementary figure panels A–C; `SFx_caption.txt`. |
