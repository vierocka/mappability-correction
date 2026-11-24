# Excluded MANE transcripts

117 per sample (115 in JY) leave quantification and TPM normalisation.
Applied by `scripts/apply_uniqueness_and_average.py` (rule 6). Tables list every ID.

| file | n | reason |
|---|---|---|
| `1_non_coding.tsv` | 70 | no protein → no peptides (RNU snRNAs, lncRNAs, RMRP …) |
| `2_pseudoautosomal_PAR.tsv` | 17 | identical X/Y copies; funique = 0; one value cannot be split between copies |
| `3_alt_copies_folded_into_primary.tsv` | 7 | alt-haplotype duplicate of a primary MANE gene; **counts added to the primary**, alt row removed |
| `4_population_variable_alt_haplotype.tsv` | 20 | alt-haplotype only, no primary copy; presence/copy number varies (KIR, LILRA3, GSTT1 …). HLA-DRB3/-DRB4 kept in JY (DR52/DR53) |
| `5_no_unique_primary_locus.tsv` | 3 | fix-patch model without one unique primary locus (NBPF1, PRAMEF22, FP885538.1) |
| `6_reincluded_non_primary_f1.tsv` | 34 | **kept** (f = 1): 31 fix-patch genes with one primary locus (MAPQ 60) + FUT3, CASP12, OR12D1 |

Effect on library share (raw TPM, % of sample): see `results/*_steps.tsv` and `*_composition.png`.
Largest in total-RNA HCT116 (non-coding 58 % of raw TPM from fcM counts, 40 % from unique counts).

Sources: `resources/non_primary_MANE_transcripts.tsv` (curated, `excluded`/`notes` columns),
`resources/non_coding_MANE_transcripts_no_protein.tsv`, constants `PAR_GENES`,
`FOLD_INTO_PRIMARY`, `CELL_SPECIFIC_INCLUDE` in the correction script.

Rationale: MANE Select is one transcript per gene on GRCh38 [Morales 2022]; alt haplotypes and
patches are outside the primary assembly used for funique (Ensembl primary, nf-core recommendation).
Multi-copy / identical sequences cannot be resolved by short reads [Robert & Watson 2015; Deschamps-Francoeur 2020; Derrien 2012].
