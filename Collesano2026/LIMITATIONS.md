# Limitations

What this pipeline does not do, does not measure, or measures with a known bias. Each entry gives the
number that bounds it and where that number comes from.

## The self model

- **One MANE Select protein per gene**, so abundance is effectively gene-level. All RefSeq isoforms give
  19,866 genes and 12.09 M unique 9-mers against 19,289 and 10.94 M for MANE (+10 %): extra isoforms add
  few new 9-mers but split the expression. [LESSONS](LEARNED_LESSONS.md)
- **Non-canonical and cryptic peptides are out of scope.** Reported shares of the ligandome range from
  ~10 % (Laumont 2016) to 1.7 % in a 772-sample reanalysis (2023), and cryptic ORFs yield MHC-I peptides
  ~5× more efficiently per translation event (Ouspenskaia 2021). This bounds achievable recall.
- **Not modelled:** peptides revealed when NMD is inhibited (Vendramin 2026), peptides from translated
  pseudogenes (Cheetham 2020; Lin 2025), translation efficiency (Madern 2025), mRNA–protein discordance
  (Upadhya 2022). [LITERATURE](LITERATURE.md)
- **Only HLA-A/B/C are personalised.** The rest of the proteome is reference sequence, so the cell lines
  differ in expression, not in peptide sequence, outside HLA.

## The mappability correction

- **It is triggered by sequence similarity, not by pseudogene expression.** A pseudogene only steals reads
  if it is transcribed. Three of ten chr6 pseudogenes examined have lost the parent promoter entirely
  (HLA-DRB6, -DPA2, -DQB3) and should contribute little, yet similarity alone would trigger the
  correction. The error in either direction is bounded by 1 − UF.
  [pseudogene_problem](pseudogene_problem/README.md)
- **The floor for f = 0 is an arbitrary constant** (0.001, `ZERO_FLOOR` in
  `scripts/apply_uniqueness_and_average.py`). At most 45 protein-coding transcripts can reach it — those
  with UF = 0 *and* pseudogene evidence, the only ones eligible for the correction rule.
- **Decoys are not implemented.** Pseudogene competition is inferred from similarity and scaled afterwards,
  rather than measured by putting pseudogene sequences in the index.
  [decoys.md](pseudogene_problem/decoys.md)
- **funique is computed from exon ± 50 bp genomic windows**, so no junction-spanning reads are simulated.
  Precomputed tables exist for 75 bp and 100 bp reads; other lengths are simulated on demand.

## Quantification

- **Unstranded counting** (`featureCounts -s 0`) on libraries that are stranded, so antisense and nested
  genes in the APM locus cannot be separated. Open; it needs the BAMs.
  [PHENOMENA row 6](PHENOMENA.md)
- **Library chemistry differs between the lines and is not removed by non-coding filtering.** HCT116 is
  total RNA (Ribo-Zero), JY and HCC1937 are poly(A). After removing non-coding transcripts and rescaling,
  the 69 replication-dependent histone genes — the one protein-coding family without a poly(A) tail — still
  hold 1.70 % of coding TPM in HCT116 against 0.036–0.061 % in the poly(A) lines (H4C3: ~1,150 vs ~1 TPM).
  Z-scaling in later steps equalises magnitude but not the direction of this bias.
  [tests/library_composition](tests/library_composition/README.md)
- **The transcript TPM ≥ 1 cut-off** removes 24–34 % of 9-mers and costs 0.7–1.9 points of MS recall. No
  published RNA → immunopeptidomics method applies an expression cut-off when building a repertoire.
  [tests/filter_impact](tests/)
- **Occurrence weighting (× occurrences) is untested here**, and is near-irrelevant in practice: 99.6 % of
  peptides occur once in their protein (pepX).

## Validation

- **The correction is not validated against mass spectrometry.** Global metrics cannot test it, because its
  effect is confined to the extreme tail (31–36 % of the top 0.1 % of net sums): over 9.4 M 9-mers it moves
  |ΔAUC| ≤ 0.0004, recall +0.00 to +0.04 points and Spearman ≥ 0.9977. The tail-restricted test has been run
  and is also null — PPV and recall inside the top 0.1–1 % change by one to three peptides with inconsistent
  sign, and against abundance-matched non-ligands the rank shift of MS ligands is indistinguishable
  (p = 0.76–0.81). A gene-level test — do demoted transcripts carry fewer ligands than their uncorrected
  abundance implies — is suggestive but not established: pooled rate ratio 0.75 (0.57–0.97) under a negative
  binomial fit and 0.85 (0.67–1.08) under a model-free stratified estimator, with HCC1937 pointing the other
  way. The limit is structural: MS ligands sit at a median 85th–87th percentile of net sum, only 1–2 % reach
  the top 0.1 %, and just 7–31 ligands per line land on demoted transcripts.
  [tests/tail_validation](tests/tail_validation/README.md)
- **Recall is a sanity check, not a discriminating metric.** It sits at 97.7–98.1 % (98.6–98.9 % without the
  TPM ≥ 1 cut-off), because nearly every MS ligand comes from a gene detected in these bulk RNA-seq libraries.
  The × funique step moves it by +0.00 to +0.04 points. Non-coding exclusion is the one step that moves it
  materially, and only in the Ribo-Zero line: HCT116 96.64 → 97.68. [tests/filter_impact](tests/)
- **No held-out data.** Filters and thresholds were chosen while looking at the Bassani-Sternberg 2015
  outputs, and the same study is then used as validation. The per-line MS tables contain only ligands already
  present in the source repertoire (100 % recall by construction), so they cannot test recall.
- **The filter comparison changes its own evaluation set**, since filtered peptides also leave the negative
  set, so AUC differences between cut-off variants are partly mechanical.

## Data and scope

- **Three cell lines**: two carcinoma-derived (HCC1937, breast; HCT116, colorectal) and JY, an
  EBV-transformed B-lymphoblastoid line. HCT116 has lane replicates only, not biological replicates.
- **HLA typed from short reads only** — two typers plus Cellosaurus reconciliation, no orthogonal typing.
- **Not compared with long reads** (isoform truth) or with proteomics.
- **Degradation is unmeasured, not excluded.** RIN is unknown for public SRA runs and 3′ bias can mimic a
  shorter isoform; with one transcript per gene the effect reduces to gene level. [LESSONS](LEARNED_LESSONS.md)
- **The server workflow has not been executed end to end here.** A dry run verifies that the DAG connects
  and that inputs exist; it does not execute the rules. [HOW2RUN](HOW2RUN.md)
