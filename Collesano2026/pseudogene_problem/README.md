# The pseudogene problem

Why the uniqueness factor (funique) exists, and why it is applied to one cause of low uniqueness and not the
others. Reads from a pseudogene locus that is absent from the MANE reference have nowhere else to go: they are
assigned to the functional parent gene, which inflates its count. The correction `count × funique` removes that
share. This folder holds the evidence for the problem and the limits of the correction.

## Low uniqueness has four causes, and only one of them inflates

Of the 1,444 protein-coding MANE transcripts with low uniqueness
(`../resources/funique/low_uniqueness_MANE_transcripts.tsv`; non-coding excluded as in the pipeline):

| cause | transcripts | median UF (L100) | effect on the count | corrected? |
|---|---|---|---|---|
| **pseudogenes** | 737 (329 with no other cause) | 0.71 | **inflation** — reads from loci outside MANE land on the parent | **yes, × funique** |
| paralog in MANE | 768 | 0.61 | deflation — both copies are in the reference and compete | no |
| duplicate in MANE | 461 | 0.10 | deflation — identical copies split or zero the signal | no; excluded or folded |
| self-homology | 52 | 0.55 | deflation — internal tandem repeats (NBPF, MUC, EPPK1) | no |
| none detected | 101 | 0.84 | — | no |

The distinction matters: multiplying by funique is only valid when the missing competitor is *outside* the
reference. When the competitor is inside MANE, the reads are already divided between the two records, and
multiplying again would correct twice. This is why the correction rule also requires
`shared = 1 − fc/fcM < 0.10` (no competition inside MANE) before it fires.

## Evidence 1: expression sits on the pseudogene-bearing transcripts

`inflation_evidence.py` ranks the flagged transcripts by observed count. Across all five libraries, almost every
one of the most abundant is pseudogene-bearing, against a 51 % base rate:

| library | top 15 with pseudogene evidence | binomial p | median count (all) | median count (top 15) | the remainder |
|---|---|---|---|---|---|
| JY | 15/15 | 4 × 10⁻⁵ | 45 | 170,869 | — |
| HCC1937 rep1 | 14/15 | 6 × 10⁻⁴ | 26 | 35,566 | UBC |
| HCC1937 rep2 | 14/15 | 6 × 10⁻⁴ | 22 | 35,942 | UBC |
| HCT116 rep1 | 14/15 | 6 × 10⁻⁴ | 22 | 49,566 | UBC |
| HCT116 rep2 | 14/15 | 6 × 10⁻⁴ | 22 | 50,536 | UBC |

The one exception is the same gene every time: **UBC**, whose low uniqueness comes from its own tandem
polyubiquitin repeat (self-homology, UF 0.54–0.62), not from pseudogenes. The classification separates the two
mechanisms correctly.

Recurring across libraries: EEF1A1 (10 competing pseudogene loci, 50 named pseudogenes; 2.09 M reads in JY),
ribosomal proteins (RPL7 with 29 competing loci, RPS3A, RPLP0, RPL41, RPL26, RPL9), TMSB4X, EEF1G, PTMA, NPM1,
HNRNPA1, PPIA — and **HLA-A**.

This is why the correction matters for a net-sum repertoire even though it leaves global metrics unchanged: a
9-mer's net sum is dominated by the most abundant transcripts encoding it, and those are the affected ones. The
share of 9-mers corrected rises from 1.2–1.7 % in the bottom half of net sums to 31–36 % in the top 0.1 %
(`../tests/inflation_tail/`).

**Non-coding transcripts must be excluded first.** In the total-RNA HCT116 libraries the multicopy snRNAs
(RNU1-1…4, RNU2-1/2; UF 0–0.83) otherwise take six of the top seven places and the enrichment disappears
(8/15, p = 0.53). They are a deflation artefact of multicopy annotation, not inflation, and carry no peptides.
`inflation_evidence.py` removes them by default; `--keep-non-coding` reproduces the effect.

```bash
python3 inflation_evidence.py -n 15 --run fcM_JY_SRR364065      --uf UF_L100
python3 inflation_evidence.py -n 15 --run fcM_HCC1937_SRR31477874 --uf UF_L75
```

## Evidence 2: not every pseudogene competes

A pseudogene only steals reads if it is transcribed. [`tss_utr_degeneration/`](tss_utr_degeneration/README.md)
compares ten chr6 MHC pseudogenes with their functional parent, segment by segment (GENCODE v49, GRCh38):

- **HLA-H, -J, -K** keep the parent's promoter, CDS and 3′UTR at 85–94 % — near-complete copies that can be
  transcribed, and are derepressed in cancer ([HLApseudo_MAFFT](../../HLApseudo_MAFFT/README.md)).
- **HLA-L, -W, -DPB2, MICE** keep the CDS (53–84 %) but have lost much of the 3′UTR (33–66 %).
- **HLA-DRB6, -DPA2, -DQB3** retain no trace of the parent promoter or 5′UTR; DPA2 and DQB3 also lack the 3′UTR.

So sequence similarity alone overstates the competition. This is the main known weakness of the correction: it
is triggered by similarity, while the inflation it corrects requires the pseudogene to be *expressed*. See
[LIMITATIONS](../LIMITATIONS.md).

## What is not implemented

[`decoys.md`](decoys.md) sets out the alternative: add pseudogene sequences to the index as decoys and let
fractional assignment separate the reads, instead of scaling afterwards. That measures competition directly
rather than inferring it from similarity, and it is the principled fix for the weakness above. Not run.

## Files

| file | content |
|---|---|
| `inflation_evidence.py` | ranks flagged transcripts by count; binomial enrichment against the base rate |
| `pseudogene_inflation_top_*.tsv` | top 15 per library, with UF, competing loci and named pseudogenes |
| `tss_utr_degeneration/` | chr6 pseudogene vs parent, by promoter / 5′UTR / CDS / 3′UTR |
| `decoys.md` | design for the decoy-index test (not run) |

Related: [`../APM_locus/`](../APM_locus/README.md) (chr6 geometry), [`../PHENOMENA.md`](../PHENOMENA.md) row 3,
[`../tests/inflation_tail/`](../tests/README.md), [`../LITERATURE.md`](../LITERATURE.md) (Cheetham 2020, Lin 2025).
