# Library chemistry: does removing non-coding RNA make the lines comparable?

JY and HCC1937 are poly(A)-selected; HCT116 is total RNA (Ribo-Zero). TPM is compositional, so a library
whose mass sits in non-coding RNA places every protein-coding gene on a different scale. The pipeline removes
70 non-coding MANE transcripts and renormalises. This test asks whether that is sufficient.

**It is not.** The rescaling fixes the overall scale but leaves a protein-coding family off by three orders of
magnitude.

```bash
python3 compare_polya_ribozero.py --counts-root <data root> --gene-map mane_gene_map.tsv
```

## Composition before and after rescaling

| library | chemistry | non-coding % of raw TPM | rescale factor | top-10 % after | top-100 % after | median TPM | % above TPM 1 | RD histones, % of coding TPM |
|---|---|---|---|---|---|---|---|---|
| JY | poly(A) | 0.17 | 1.002 | 11.1 | 41.1 | 11.0 | 58.9 | **0.036** |
| HCC1937 r1 | poly(A) | 0.07 | 1.001 | 6.2 | 25.1 | 14.0 | 65.4 | **0.061** |
| HCC1937 r2 | poly(A) | 0.07 | 1.001 | 6.3 | 25.4 | 13.7 | 64.8 | **0.056** |
| HCT116 r1 | Ribo-Zero | 40.43 | 1.679 | 12.7 | 43.7 | 10.0 | 63.5 | **1.707** |
| HCT116 r2 | Ribo-Zero | 40.45 | 1.679 | 12.7 | 43.7 | 10.0 | 63.6 | **1.701** |

(Unique counts, `fc`. The 40 % non-coding share matches [PHENOMENA](../../PHENOMENA.md) row 5, which reports
58 % for `fcM` and 40 % for unique counts.)

After rescaling the summary statistics do converge: median TPM 10–14, share above TPM 1 within 59–65 %,
top-100 concentration 25–44 %. On those measures the libraries look comparable.

## The chemistry marker that survives

Replication-dependent histones are the one protein-coding family with **no poly(A) tail**. Poly(A) selection
depletes them; Ribo-Zero retains them. They are therefore a marker of library chemistry that is independent of
which genes a cell line happens to express. (Replication-*independent* variants — H1-0, H1-10, H2AZ, H2AX,
H3-3, H3-5, macroH2A — are polyadenylated and are excluded from the set.)

After non-coding removal and rescaling, 69 replication-dependent histone genes hold **1.70 % of protein-coding
TPM in HCT116 against 0.036–0.061 % in the poly(A) lines: a 28–47-fold difference.** Per gene it is larger:

| gene | HCC1937 r1 | HCC1937 r2 | JY | HCT116 r1 | HCT116 r2 |
|---|---|---|---|---|---|
| H4C3 | 1.0 | 0.5 | 0.7 | 1,149 | 1,146 |
| H1-4 | 1.9 | 1.3 | 0.7 | 1,091 | 1,105 |
| H2AC20 | 1.2 | 1.7 | 0.6 | 904 | 879 |
| H4C2 | 1.3 | 0.3 | 0.2 | 713 | 702 |
| H3C10 | 6.6 | 5.0 | 3.3 | 657 | 675 |
| H3C2 | 2.1 | 1.0 | 2.0 | 575 | 599 |
| H1-2 | 85.6 | 79.9 | 111.1 | 1,607 | 1,558 |

H4C3 differs by ~1,500-fold between HCT116 and the poly(A) lines. This is not biology: HCT116 is not making a
thousand times more histone H4 than JY. It is the poly(A) step.

## Consequences

1. Non-coding removal and rescaling do not remove the chemistry difference. It is concentrated in a family that
   encodes peptides and therefore enters the repertoire. Z-scaling across samples, applied in later steps, can
   equalise magnitude but not the direction of a bias caused by different library preparation chemistries.
2. **Histone-derived 9-mers are a systematic artefact in HCT116.** Histones are short (~100–220 aa), so
   69 genes contribute on the order of 10⁴ 9-mers out of 9.4 M; carrying ~1.7 % of coding TPM, those peptides
   land in the extreme tail of the net-sum distribution — the same region where the mappability correction acts
   and where the repertoire is most sensitive.
3. Histones are independently among the hardest genes to quantify: intronless, nearly identical within
   clusters, no poly(A) 3′UTR to disambiguate, and among the ~40 genes that UTRs fail to rescue
   ([LESSONS §3](../../LEARNED_LESSONS.md)). Chemistry bias and mappability bias coincide on the same family.

**Recommended handling:** restrict quantitative comparisons to within a cell line, and either exclude
replication-dependent histones from cross-line statements or report them separately. The current pipeline does
neither; see [LIMITATIONS](../../LIMITATIONS.md) §3.

**Caveat.** The three lines differ in both chemistry and biology, so a matched poly(A)/Ribo-Zero pair from one
line would be needed to separate the two cleanly. The histone family is used here precisely because its
polyadenylation status makes the chemistry interpretation unambiguous, but the magnitude for any individual
gene still mixes the two.

## Files

| file | content |
|---|---|
| `compare_polya_ribozero.py` | the analysis |
| `library_composition.tsv` | the table above |
| `histone_tpm_by_library.tsv` | all 69 replication-dependent histones, TPM after rescaling |
| `mane_gene_map.tsv` | transcript → gene symbol, from the MANE FASTA headers |

## Filter-step tables

`filter_steps_*.tsv` (from `../../scripts/report_filter_steps.py`, five per-sample inputs) back
[PHENOMENA row 5](../../PHENOMENA.md): non-coding is 57.8 % of raw fcM TPM in HCT116 against
0.08–0.18 % in the poly(A) lines, and the total rescale from all exclusions plus × funique is
×1.0371 (JY), ×1.0305 (HCC1937) and ×2.4448 (HCT116). Those rescale factors are reproduced
independently by [tests/tail_validation](../tail_validation/README.md).

| file | content |
|---|---|
| `filter_steps_steps.tsv` | every exclusion and correction step, per sample: transcripts, counts, % of raw counts and TPM |
| `filter_steps_composition.tsv` | biotype composition before and after |
| `filter_steps_shift.tsv` | median log2 TPM shift per sample (the rescale factor) |
| `filter_steps_groups.tsv` | transcripts removed per exclusion rule |
| `filter_steps_top_excluded.tsv` | the transcripts carrying the most excluded TPM |
| `filter_steps_steps.png`, `filter_steps_composition.png` | the same as figures |
