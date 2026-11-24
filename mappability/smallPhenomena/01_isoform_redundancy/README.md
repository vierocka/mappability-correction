# 01 — Isoform redundancy  [v2]

**Central question:** When a gene has multiple annotated CDS isoforms, or belongs to a
paralog family, how much signal can short-read RNA-seq tools actually resolve — and does
the UTR-CDS-UTR transcript structure help?

---

## Biological context

The deduplicated RefSeq CDS FASTA contains all canonical coding sequences after exact-
duplicate removal (~93,000 entries across ~19,000 genes). Many genes have multiple
isoforms — alternative splicing produces CDS variants that share most exons. When a
31-mer (k=50 used here; kallisto pseudoalignment default is k=31) appears in two isoforms of the same gene, reads from
either isoform land in the same equivalence class; EM redistribution splits the signal
— it does not distinguish them.

The CDS-to-CDS self-map (STAR simulated reads back to the same FASTA) quantifies this
directly: a uniqueness factor (UF) near zero means essentially all reads from that
sequence also match another sequence in the reference.

**UTR rescue:** switching from a CDS-only index to a full UTR-CDS-UTR MANE index rescues
~90% of genes (CDS frac_unique ≈ 3% → MANE RNA frac_unique ≈ 94%). The remaining ~10%
(~1,922 genes) are unrescued. Why_UTRs_fail_to_rescue_nearly2000annotGenes.txt explains
the five mechanistic categories.

---

## Analyses

### A. Isoform table (`isoform_table.tsv`)
One row per gene from the dedup CDS FASTA. Gene name, isoform count, lengths, protein
IDs. Sequence-independent; run once.

### B. Functional group enrichment (`isoform_groups.tsv`)
Joins isoform counts with:
- **Pseudogene count** — GENCODE v49 gene entries with pseudogene biotypes, matched to
  coding parent by gene_name prefix (see note below)
- **Functional group membership**: kinase, TF, cell cycle, extracellular, functional
  blacklist (from `~/Reference/MANE/*_enst.txt`)
- **ID bridge**: `../results/id_conversion_table.tsv` maps NP_ → ENST → gene_symbol

Questions answered: do multi-isoform genes accumulate more pseudogenes? Are kinases /
TFs enriched among single-isoform genes?

#### Note: how pseudogene–gene affiliation is defined here

**What we look for in GENCODE v49:**
Any feature-level entry (`gene` record) carrying a pseudogene `gene_type`, specifically:
`processed_pseudogene`, `unprocessed_pseudogene`, `transcribed_processed_pseudogene`,
`transcribed_unprocessed_pseudogene`, `transcribed_unitary_pseudogene`,
`unitary_pseudogene`, `rRNA_pseudogene`, `IG_V_pseudogene`, `TR_V_pseudogene`,
`TR_J_pseudogene`.

**How we define affiliation (parent→pseudogene):**
GENCODE does not always populate the `parent_gene_id` attribute for pseudogenes (it is
present for retropseudogenes derived from a specific source transcript, but absent for
many duplicated pseudogenes). We therefore use a **gene_name prefix heuristic**:

A pseudogene is affiliated to coding gene G if:
1. The pseudogene `gene_name` starts with the same uppercase string as G, AND
2. The next character after G's name is `P`, `_`, or `-`
   (e.g. `OR4F5` → `OR4F5P1`, `TP53` → `TP53BP1` would NOT match because `B` ≠ `P/_/-`)

This catches the canonical HGNC pseudogene naming convention (gene + `P` + integer, or
gene + `-ps` style) while avoiding false positives between genes sharing a short prefix
(e.g. `ZNF` would not absorb `ZNFX1`). The longer coding gene name is tested first to
prevent `GP1BA` from capturing `GP1BB` as a pseudogene.

**What this count means and does NOT mean:**
- `pseudogene_count > 0` means at least one GENCODE-annotated pseudogene has a name
  that looks like it derived from this coding gene. It is a reasonable proxy for "how
  many copies were duplicated and then became non-functional."
- It is **not** a curated homology call. A pseudogene named `OR4F5P1` is affiliated to
  `OR4F5` by HGNC convention, but the actual sequence donor may be a different OR family
  member — GENCODE `parent_gene_id`, when present, is the authoritative source.
- Some pseudogenes have no obvious parent by name (e.g. unnamed processed pseudogenes);
  these are not counted here.
- Processed pseudogenes (retropseudogenes) often have very high CDS identity to their
  parent (they derive from a reverse-transcribed mRNA); duplicated pseudogenes have lower
  identity and may carry their own UTR divergence. The `pseudogene_count` column does not
  distinguish these two classes — use the GENCODE `gene_type` field for that.

### C. 31-mer unique fraction and unique stretch (`isoform_kmers_cds.tsv`, `isoform_kmers_mane.tsv`)
For each isoform/transcript: fraction of 31-mers unique to that sequence alone,
plus the length of the longest contiguous run of unique 31-mer positions.

**k-mer uniqueness (frac_unique):** fraction of all 31-mers in the sequence not
present in any other sequence in the same reference set.

**Unique stretch:** the longest run of consecutive positions i where seq[i:i+31]
is a unique 31-mer. Each position in the run adds 1 bp; a run of N positions spans
N + 30 bp. One unique 31-mer (N=1) spans 31 bp — the minimum reportable case.
Target: N ≥ 20 consecutive positions → stretch ≥ 50 bp. A 50 bp unique region
tolerates short indels or mismatches that shift k-mer positions by a few bases
while still guaranteeing ≥ 1 unique 31-mer within any overlapping read.

New columns in output files: `max_unique_stretch_bp`, `max_stretch_positions`,
`stretch_start_nt`, `stretch_reaches_50bp`.

Computed twice:
- **CDS-only** (`isoform_kmers_cds.tsv`): within-gene uniqueness, dedup CDS sequences.
  The floor: what is distinguishable from coding sequence alone.
- **MANE RNA** (`isoform_kmers_mane.tsv`): cross-gene uniqueness over all MANE
  transcripts (5'UTR + CDS + 3'UTR). frac_unique < 0.90 defines the UTR-fail
  genes used as input to Sections F and G.

The delta between the two reveals how much UTRs and splice-junction k-mers add.

### D. Per-run STAR UF + kallisto (22 conditions × 2 indices + v8 proper)

One subfolder per simulation variant under `runs/`. Each contains:
- `vs_dedup_cds/isoform_mappability.tsv` — STAR UF + kallisto on the NP_ CDS index
- `vs_dedup_cds/reads_by_nisoform.tsv`   — same, binned by n_isoforms
- `vs_MANE/isoform_mappability.tsv`       — STAR UF + kallisto on the MANE index
- `vs_MANE/reads_by_nisoform.tsv`         — same, binned by n_isoforms

v8 proper-pair runs have their own layout (ENST-keyed, no isoform split):
- `v8_proper_*/gene_uf.tsv`              — gene-level STAR UF, proper-pair geometry

#### Read conditions

| Label | Mode | Read len | Read source | STAR UF | Kallisto |
|---|---|---|---|---|---|
| `dedup_cds_reads_PE_L75bp`      | PE | 75 bp  | dedup CDS only | — | ✓ |
| `dedup_cds_reads_PE_L100bp`     | PE | 100 bp | dedup CDS only | ✓ | ✓ |
| `dedup_cds_reads_PE_L150bp`     | PE | 150 bp | dedup CDS only | — | ✓ |
| `dedup_cds_reads_PE_L200bp`     | PE | 200 bp | dedup CDS only | — | ✓ |
| `dedup_cds_reads_SE_L75bp`      | SE | 75 bp  | dedup CDS only | — | ✓ |
| `dedup_cds_reads_SE_L100bp`     | SE | 100 bp | dedup CDS only | ✓ | ✓ |
| `dedup_cds_reads_SE_L150bp`     | SE | 150 bp | dedup CDS only | — | ✓ |
| `dedup_cds_reads_SE_L200bp`     | SE | 200 bp | dedup CDS only | — | ✓ |
| `MANE_RNA_reads_PE_L75bp`       | PE | 75 bp  | MANE full RNA  | — | ✓ |
| `MANE_RNA_reads_PE_L100bp`      | PE | 100 bp | MANE full RNA  | — | ✓ |
| `MANE_RNA_reads_PE_L150bp`      | PE | 150 bp | MANE full RNA  | — | ✓ |
| `MANE_RNA_reads_PE_L200bp`      | PE | 200 bp | MANE full RNA  | — | ✓ |
| `MANE_RNA_reads_SE_L75bp`       | SE | 75 bp  | MANE full RNA  | — | ✓ |
| `MANE_RNA_reads_SE_L100bp`      | SE | 100 bp | MANE full RNA  | — | ✓ |
| `MANE_RNA_reads_SE_L150bp`      | SE | 150 bp | MANE full RNA  | — | ✓ |
| `MANE_RNA_reads_SE_L200bp`      | SE | 200 bp | MANE full RNA  | — | ✓ |
| `MANE_exonflank_reads_PE_L75bp` | PE | 75 bp  | exon±50 bp flank | — | ✓ |
| `MANE_exonflank_reads_PE_L100bp`| PE | 100 bp | exon±50 bp flank | — | ✓ |
| `MANE_exonflank_reads_PE_L150bp`| PE | 150 bp | exon±50 bp flank | — | ✓ |
| `MANE_exonflank_reads_SE_L75bp` | SE | 75 bp  | exon±50 bp flank | — | ✓ |
| `MANE_exonflank_reads_SE_L100bp`| SE | 100 bp | exon±50 bp flank | — | ✓ |
| `MANE_exonflank_reads_SE_L150bp`| SE | 150 bp | exon±50 bp flank | — | ✓ |

**v8 proper-pair (150 bp PE, real 300 bp fragment geometry):**

| Label | Reference |
|---|---|
| `v8_proper_genomic`   | Full GRCh38.p14 genome |
| `v8_proper_dedup_cds` | Deduplicated CDS reference |
| `v8_proper_MANE`      | MANE Select transcriptome |

**Index comparison:**
Each kallisto run is processed against two indices:
- `vs_dedup_cds` — NP_ CDS targets only (no UTRs); reveals isoform-level CDS ambiguity
- `vs_MANE`      — ENST targets, full UTR-CDS-UTR MANE transcripts; reveals how much
                   UTR k-mers reduce multi-mapping

**Key finding (from earlier dedup_cds_reads runs):**
Read length has minimal effect on kallisto entropy: 0.9957 (75 bp) → 0.9812 (200 bp)
for 2-isoform genes, < 0.002 difference. PE = SE at every length.
The bottleneck is the index k-mer topology (kallisto default k=31), not read length.

### E. MAFFT pairwise CDS isoform homology (`isoform_homology.tsv`)
Pairwise % sequence identity between isoforms of the same gene. Multi-isoform genes
with all isoforms ≤ 4,000 nt, capped at 500 genes for runtime. Run once.

### F. UTR-fail equal-signal statistical tests (`utr_fail_stats.tsv`)
For each of the ~1,922 UTR-fail genes (frac_unique < 0.90 in isoform_kmers_mane.tsv),
tests whether the kallisto signal is uniformly distributed across CDS isoforms.

Tests per gene × condition (across all vs_dedup_cds runs):
- **χ² goodness-of-fit vs uniform** (scipy.stats.chisquare on est_counts)
- **Normalised Shannon entropy** H/log₂(n_isoforms) ∈ [0, 1]
- **Coefficient of variation** (std/mean of est_counts across isoforms)

Aggregated across all 22 conditions; verdict column:
`uniform_signal` | `mostly_uniform` | `partially_resolved` | `no_data`

Genes with `uniform_signal` cannot be isoform-resolved by any short-read approach.

**Note:** ANOVA per functional group is a separate follow-up step:
```r
df <- read_tsv("utr_fail_stats.tsv")
aov(mean_entropy ~ functional_groups, data = df) |> summary()
```

### G. Paralog MAFFT (`paralog_homology_mane.tsv`, `paralog_homology_isoforms.tsv`)
MAFFT pairwise distance between UTR-fail genes in the same paralog family.

Families defined by HGNC symbol prefix (e.g. OR, ZNF, PCDHA, HIST1H2A).
Size filter: 2 ≤ n_members ≤ 99.

- **Option A** (`paralog_homology_mane.tsv`): one MANE canonical sequence per gene
  (full UTR-CDS-UTR from MANE RNA FASTA). One row per gene pair.
- **Option B** (`paralog_homology_isoforms.tsv`): all GENCODE v49 isoforms per member,
  if ≤ 20 isoforms per family member. Captures isoform × isoform cross-gene distances.

ID cross-verification across all four sources written to `id_verification.tsv`.

### H. Actionable gene lists
Derived from Sections C, F, and G.

- **`gene_level_recommendation.tsv`**: genes where isoform redundancy makes CDS
  isoform signal indistinguishable (uniform entropy + UTR rescue < 0.90).
  Recommendation: collapse all CDS isoforms to a single gene-level TPM sum.

- **`cluster_level_recommendation.tsv`**: genes in paralog families where
  max pairwise MAFFT identity ≥ 85% (Option A) AND UTR rescue insufficient.
  Recommendation: aggregate all family members into a paralog-cluster TPM sum;
  do not interpret individual gene TPMs.

---

## Output files

| File | Section | Produced | Description |
|---|---|---|---|
| `isoform_table.tsv` | A | once | gene → n_isoforms, lengths, protein IDs |
| `isoform_groups.tsv` | B | once | + pseudogene count + functional group membership |
| `isoform_kmers_cds.tsv` | C | once | per-isoform unique 31-mer fraction + unique stretch (CDS) |
| `isoform_kmers_mane.tsv` | C | once | per-transcript unique 31-mer fraction + unique stretch (MANE RNA) |
| `isoform_homology.tsv` | E | once | pairwise MAFFT identity, CDS isoforms, same gene |
| `runs/{label}/vs_dedup_cds/isoform_mappability.tsv` | D | per run | STAR UF + kallisto vs CDS index |
| `runs/{label}/vs_dedup_cds/reads_by_nisoform.tsv` | D | per run | binned by n_isoforms |
| `runs/{label}/vs_MANE/isoform_mappability.tsv` | D | per run | STAR UF + kallisto vs MANE index |
| `runs/{label}/vs_MANE/reads_by_nisoform.tsv` | D | per run | binned by n_isoforms |
| `runs/v8_proper_*/gene_uf.tsv` | D | per run | v8 proper-pair STAR UF (ENST-keyed) |
| `utr_fail_stats.tsv` | F | once | χ², entropy, CV for ~1,922 UTR-fail genes |
| `paralog_homology_mane.tsv` | G | once | paralog pairwise identity (MANE canonical) |
| `paralog_homology_isoforms.tsv` | G | once | paralog pairwise identity (GENCODE all isoforms) |
| `id_verification.tsv` | G | once | symbol consistency across CDS/MANE/GENCODE/conv_table |
| `gene_level_recommendation.tsv` | H | once | genes to collapse to gene level |
| `cluster_level_recommendation.tsv` | H | once | genes to aggregate at cluster level |

---

## Column glossary

| Column | Meaning |
|---|---|
| `frac_unique` | Fraction of 31-mers unique to this sequence (not shared with any other in reference) |
| `max_unique_stretch_bp` | Length in bp of the longest contiguous run of unique 31-mer positions (N positions → N+30 bp) |
| `max_stretch_positions` | Number of consecutive positions in the longest unique run |
| `stretch_start_nt` | 0-based start position of the longest unique stretch in the transcript |
| `stretch_reaches_50bp` | 1 if max_unique_stretch_bp ≥ 50 (≥ 20 consecutive unique 31-mer positions), else 0 |
| `star_uf` / `star_uf_mean` | STAR uniqueness factor: fraction of simulated reads mapping back with MAPQ=255 |
| `kallisto_tpm_entropy` | Normalised Shannon entropy of TPM distribution across CDS isoforms [0=peaked, 1=uniform] |
| `kallisto_tpm_cv` | Coefficient of variation of TPM across CDS isoforms |
| `kallisto_est_counts_sum` | Sum of EM-estimated fractional read counts across isoforms |
| `n_with_kallisto` | Number of CDS isoforms receiving ≥ 1 kallisto count |
| `mean_entropy` (Section F) | Mean normalised entropy across all 22 run conditions |
| `frac_chi2_nonsig` | Fraction of conditions where χ² cannot reject uniform signal (p > 0.05) |
| `signal_distribution` | Verdict: uniform_signal / mostly_uniform / partially_resolved / no_data |
| `max_paralog_identity` | Highest pairwise % identity with any family member (MAFFT Option A) |
| `recommendation` | collapse_to_gene_level OR aggregate_at_cluster_level |

See `../README.md` for the full project-level glossary.

---

## Run history — `runs/batch1_multimapmax40/`

All pre-2026-07-08 STAR runs under `runs/` (every subfolder listed in the table
above) were computed with `--outFilterMultimapNmax 40`. They have been moved
into `runs/batch1_multimapmax40/` to preserve them, rather than being
overwritten, since the whole mappability pipeline switched to
`--outFilterMultimapNmax 20` (matching nf-core/rnaseq's ENCODE3-standard
`star_rsem` preset — see `../../Ref/recommendations_nf.txt`) on 2026-07-08.

**Why this matters here specifically:** `outFilterMultimapNmax` caps how many
loci STAR will consider before giving up and calling a read unmapped/ambiguous
rather than reporting any alignment. For isoform-redundant genes and paralog
families — exactly what this analysis studies — a lower cap (20 vs. 40) could
plausibly change `star_uf` for any transcript with more than ~20 competing
loci, independent of the isoform-redundancy/k-mer effects this analysis is
actually trying to isolate. Any new `runs/` output going forward is on the
`--outFilterMultimapNmax 20` setting — **compare against
`batch1_multimapmax40/` once both exist to check whether this threshold
change had a material impact** on `star_uf`, `isoform_mappability.tsv`, or the
downstream recommendation files (Sections D, H) before treating results
across the two batches as comparable.

## Running

```bash
cd <SELF_NONSELF>/Reference/mappability-correction/mappability/smallPhenomena/01_isoform_redundancy
python3 analyze.py
```

Requirements: `scipy` (for χ² tests in Section F), `mafft` in PATH (for Sections E and G).
To skip MAFFT, the script continues; only Sections E and G output files will be absent.
```bash
pip install scipy
# mafft: apt install mafft  OR  conda install -c bioconda mafft
```

Sections A–E run in ~20–30 min (dominated by the 31-mer counting pass over MANE transcripts;
slightly slower than k=31 due to larger k-mer set per sequence).
Section F: ~5 min.
Section G (MAFFT): runtime scales with family count × sequences; estimated 30–90 min
for all UTR-fail paralog families.
