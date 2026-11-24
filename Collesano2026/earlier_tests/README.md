# Earlier tests (copied)

Scripts, logs and small result tables from the exploratory runs that shaped the pipeline. They were run in the
project folders `human/RNAseq` and `human/References`, which are not part of this repository.
Local paths, cluster user, group and e-mail were replaced by `<placeholders>`. The scripts are kept as run, not
refactored. Interpretation: [../LEARNED_LESSONS.md](../LEARNED_LESSONS.md).

| folder | file | what |
|---|---|---|
| `RNAseq_isoforms/` | `hpc_kallisto.sh` | kallisto 0.50.1 quant (RefSeq GRCh38.p14 CDS index, 100 bootstraps) for all libraries |
| | `HCC1937_mapping_overview.txt` | download, kallisto index/quant and featureCounts unique vs multi summary |
| | `HCC1937_STAR_Log.final.out`, `HCC1937_featureCounts.summary` | STAR 2.7.11b on full GRCh38.p14; featureCounts 2.0.6 |
| | `dedup_kallisto_accessions.py`, `check_TAP_related_transcripts.sh` | collapse patch/alt duplicate accessions (max); TAP-related transcripts |
| | `gene_isoform_sum_max.py`, `max_vs_sum_table_6x.sh`, `max_vs_sum.R`, `plot_makie_sum_vs_max_kallisto.jl` | gene-level sum vs max isoform; fc vs kallisto scatter |
| | `JY_HCC1937_9merRepertoires.R` | isoforms per gene; two-isoform dominance (0.50 spike) |
| | `canonicalFCset_leadingIsoformKallisto_plots.R`, `HCC1937_unique9mers_per_gene_canonical_vs_leading_table.tsv` | unique 9-mers: canonical (fc) vs leading isoform (kallisto) |
| | `fc_vs_kallisto_max_merged.tsv`, `isoform_ties.py` | merged gene table; tie and agreement statistics (added 2026-11-22) |
| | `cancer_cell_lines.csv`, `gene_classes_notes.txt` | library metadata; gene classes (housekeeping, HCC1937 variants) |
| `References_peptides/` | `peptides_per_gene_MANE_vs_GRCh38p14.py`, `UniquePeptides_per_gene_MAN_vs_GRCh38p14_notes.txt`, `MANE_vs_GCF_UniqPeptides_per_gene_hist_improvedVisualisation.py` | unique 9-mers per gene: MANE vs all RefSeq isoforms (Suppl. Fig. 1B) |
| | `GRCh38p14_deduplication_note.txt` | `seqkit rmdup` of RefSeq CDS and proteins |
