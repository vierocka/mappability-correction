# Validation tables

Outputs of a workflow run on GENCODE v49 / GRCh38 with the project's existing funique table, kept
so the numbers in [../README.md](../README.md) can be checked without rerunning. Gzipped TSVs.

| file | rows | what it backs |
|---|---|---|
| `gene_classes.tsv.gz` | 78,691 genes | the per-cause counts and mean funique table; `cause` is the top-priority assignment, `all_causes` lists everything detected |
| `clusters_by_gene.tsv.gz` | per gene | 236,411 transcripts in 3'-end clusters |
| `selfhom.tsv.gz` | 254 genes | the self-homology set (NBPF20/19/10, EPPK1, CR1, LPA, NEB) |

Read directly: `pandas.read_csv("gene_classes.tsv.gz", sep="\t")`.
