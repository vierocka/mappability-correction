# Pseudogene decoys in the MANE index: why it is tricky

**Idea.** Add HLA pseudogene sequences (HLA-H, J, K, L, N, P, T, U, V, W; DRB2/6/7/8/9 …) to the
personal-MANE index. Reads that come from a transcribed pseudogene then land on the decoy and are
discarded, instead of inflating HLA-A.
**Status: not implemented.** Below are the problems a correct implementation must solve, and the test
we propose.

## Why it is not a simple fix

| # | problem | consequence | possible handling |
|---|---|---|---|
| 1 | **Annotated pseudogene models are partial.** GENCODE HLA pseudogene transcripts are 500–1100 bp; the genomic loci are longer (`HLApseudo_MAFFT` caveat) | reads from unannotated pseudogene exons or retained introns still land on HLA-A | use the genomic pseudogene locus (± flank) as decoy, not the transcript model (as in "gentrome" decoys) |
| 2 | **Ties.** In exons 2–3, a read can be equally good on HLA-A and on HLA-H | STAR gives MAPQ < 255 and fcM splits 50/50: **inflation becomes deflation** | count ties for the functional gene, splitting only reads with a better decoy hit (AS_decoy > AS_gene); test both rules |
| 3 | **Reference pseudogene vs non-reference functional allele.** The decoy carries GRCh38 HLA-H; the sample's HLA-A allele may be closer to GRCh38 HLA-H than to the GRCh38 HLA-A reference | true HLA-A reads are stolen by the decoy | decoys only make sense with a **personal** functional allele (already done); then true reads match their own allele best |
| 4 | **Pseudogenes are polymorphic and haplotype-dependent.** HLA-H, J, L … have alleles; some haplotypes lack them; the DRB pseudogene set depends on the DR haplotype (DR51/52/53) | one decoy set does not fit every sample: missing decoys leave inflation, surplus decoys can steal reads | per haplotype group, or IMGT pseudogene alleles (IPD-IMGT/HLA lists HLA-H, J, K, L … alleles) |
| 5 | **funique must change with the reference.** funique was simulated against the primary genome, which already contains the pseudogenes | decoys plus ×funique would correct twice (HLA-A ×0.89 and decoy loss) | with decoys, recompute funique against the same decoy-augmented index, or drop ×funique for genes that have decoys |
| 6 | **Library size and TPM.** Decoy counts are reads that belong to no peptide-coding transcript | if left in, decoys change every TPM; if removed, like the 117 exclusions | remove before TPM; report the decoy share per sample (a de-repression read-out) |
| 7 | **Expression is sample-specific.** Pseudogene de-repression varies (cancer hypomethylation) | an expected-sharing factor (funique) cannot follow it; decoys can, but only if 1–6 hold | this is the real gain of decoys: a sample-specific correction |
| 8 | **No ground truth in short-read data** | cannot tell a stolen read from a correct one | simulation (below); long reads spanning the inactivating mutations (PacBio PureTarget/Iso-Seq); MS cannot help (pseudogenes give no peptides) |

## Proposed test (simulation, local scale)

1. Truth: reads simulated from the personal HLA-A/B/C alleles (UTR + CDS) **and** from pseudogene loci at
   ratios 0, 5, 20 and 50 % of HLA-A, at 75/100/150 bp, with 0.1–1 % sequencing error.
2. Indexes: (a) personal MANE (current); (b) + transcript-model decoys; (c) + genomic-locus decoys (± 500 bp);
   (d) + IMGT pseudogene alleles.
3. Counting: fcM; tie rule "functional wins" vs "split".
4. Score: recovered HLA-A / truth for each design. Current pipeline = (a) × funique.
   The decoy is worth it only if the error stays near 0 across **all** pseudogene ratios; funique is right only at its implicit ratio.
5. Real data: decoy share in HCT116 (total RNA), HCC1937 and JY, and HLA-A TPM with and without decoys.

Related: Robert & Watson 2015; Deschamps-Francoeur 2020; Aguiar 2019 (personalised HLA index).
The decoy approach of salmon (Srivastava et al., Genome Biol 2020) is not in the paper's reference list.
