# HLA pseudogene MAFFT alignments

## Motivation

The MHC locus (chromosome 6p21.3) is the most complex region in the human genome for
short-read RNA-seq interpretation. Several properties compound:

**Structural complexity.** Classical and non-classical HLA genes are embedded within a
~4 Mb block that also contains dozens of pseudogenes, non-coding RNAs, overlapping
regulatory elements, and segmental duplications. Gene boundaries blur: pseudogene promoters
overlap functional gene enhancers; intronic sequences of one locus are exonic in another.

**Sequence homology.** HLA pseudogenes arose by tandem duplication from the same ancestral
gene family. Exons 2 and 3 — which encode the peptide-binding groove and are the most
clinically relevant regions — retain 70–95% nucleotide identity between functional HLA
alleles and their pseudogene relatives. A 75–150 bp read placed anywhere in this conserved
block cannot be uniquely attributed to a single locus: it maps equally well to HLA-A and
HLA-H, or to HLA-DRB1 and HLA-DRB6.

**Pseudogene transcription.** The pseudogenes here are *transcribed* unprocessed
pseudogenes — they retain promoters, intron/exon structure, and expression. In normal
tissue, many are epigenetically silenced (promoter CpG methylation). In cancer, global
hypomethylation and chromatin remodelling at 6p21 derepress them, adding a sample-specific
pool of ambiguous reads that does not exist in matched normal tissue. This is not a
constant bias; it is condition-differential and cannot cancel in differential expression
analysis.

**Inactivating mutations are sparse.** Stop codons, frameshift indels, and disrupted
splice sites occur at only 1–3 positions per pseudogene. Most of the transcript — and
therefore most reads — is intact and indistinguishable from the functional allele.

**Consequence for RNA-seq.** Standard pipelines either:
- discard multimappers (STAR default + featureCounts) → uniform signal attenuation in
  every sample; relative comparisons are internally consistent but absolute quantification
  is meaningless for HLA;
- redistribute reads by EM (kallisto, Salmon, RSEM) → inflation or deflation that is
  proportional to the pseudogene transcription level in each sample; differential
  pseudogene derepression in cancer manufactures apparent HLA fold changes that do not
  reflect protein surface expression.

**Why immune evasion conclusions require orthogonal evidence.** HLA class I loss is a
genuine immune evasion mechanism. But it cannot be inferred from RNA-seq alone when the
locus is structurally and epigenetically complex in precisely the tissue being studied.
Independent evidence should include at minimum one of: IHC or flow cytometry for HLA
surface expression; WES/WGS for loss of heterozygosity at 6p21; targeted bisulfite
sequencing at HLA pseudogene promoter CpG islands (to determine whether pseudogene
derepression is contributing to the RNA-seq signal); or long-read RNA-seq, which spans
inactivating mutations and resolves transcripts unambiguously.

**Personalised reference.** HLA is the most polymorphic locus in the human genome; each
individual carries two alleles per gene, differing by up to hundreds of substitutions.
Standard GRCh38 encodes a single mosaic haplotype. Alignment to a personalised reference
(HLA typed by HLA*LA, T1K, or arcasHLA; incorporated into the STAR index as additional
sequences) is prerequisite for allele-level quantification and substantially reduces
multimapping by providing allele-specific k-mers absent from pseudogenes.

---

## Folder architecture

```
HLApseudo_MAFFT/
│
├── run_hla_mafft.py          Main pipeline (transcript extraction + MAFFT)
├── shared_fraction.py         Sliding-window shared-sequence analysis
├── explore_commands.sh        Exploratory bash commands used during design
│
├── transcript_table.tsv       Selected transcript per gene (ID, length, biotype, rule)
├── shared_fraction.tsv        % windows shared per gene × window size × threshold
├── shared_positions.tsv       Per-position shared flag at W=100 bp, 90% (for plotting)
│
├── classI_A_clade.input.fa    Unaligned: HLA-A, E, G  +  pseudogenes H, J, L, P, V
├── classI_A_clade.aln.fa      MAFFT L-INS-i alignment, FASTA
├── classI_A_clade.aln.clustal MAFFT L-INS-i alignment, Clustal (conservation track)
├── classI_A_clade.aln.log     MAFFT log
│
├── classI_BC_clade.*          HLA-B, C, F  +  pseudogenes K, N, T, U, W
├── classII_DRB.*              HLA-DRB1, 3, 4, 5  +  pseudogenes DRB2, 6, 7, 8, 9
├── classII_DQ.*               HLA-DQA1, DQB1, DQA2, DQB2  +  pseudogene DQB3
└── classII_DP.*               HLA-DPA1, DPB1  +  pseudogenes DPA2, DPB2
```

Sequences are full spliced transcripts (5'UTR + CDS + 3'UTR) from GENCODE v49.
Functional HLA: MANE Select transcript. Pseudogenes: longest annotated transcript.
Many pseudogene transcript models are partial (they retain only a subset of exons),
appearing as leading/internal gaps in the alignment. This reflects the biological reality:
pseudogene annotation from expression evidence captures only the expressed fragments.

**Important caveat.** Pseudogene transcript lengths in GENCODE (500–1100 bp) reflect
what could be annotated from expression evidence, not the full extent of genomic sequence
similarity. The actual pseudogene loci in the genome are much longer and span regions
that include intronic sequence — which is more diverged. For mappability purposes what
matters is the exonic overlap between transcribed sequence, which these alignments
capture. The shared-sequence fractions computed by `shared_fraction.py` are therefore
lower bounds: they are based on the annotated transcribed portion of the pseudogene only.
Full genomic alignment (including unannotated exons and retained introns) would reveal
additional shared blocks. Conversely, alignments are based on GRCh38 reference alleles;
specific patient alleles may differ, and some alleles may have higher or lower sequence
identity to specific pseudogenes than the reference haplotype shown here.

---

## Scripts

**`run_hla_mafft.py`**  
Reads GENCODE v49 comprehensive transcript FASTA
(`GENCODE/gencode.v49.transcripts.fa.gz`) and MANE v1.5 RNA FASTA for MANE Select
transcript IDs. For each of the five paralog groups, selects the representative
transcript per gene (MANE Select for functional HLAs; longest pseudogene-biotype
transcript for pseudogenes), writes a per-group multi-FASTA, and runs MAFFT L-INS-i
(`--localpair --maxiterate 100`), chosen because the sequences share one conserved
domain (peptide-binding groove) flanked by divergent UTRs and partially deleted
pseudogene exons — the structural profile L-INS-i handles best. Outputs FASTA and
Clustal alignment files plus a `transcript_table.tsv` recording which transcript was
selected for each gene.

**`explore_commands.sh`**  
Bash one-liners run during design to survey GENCODE v49 HLA content, enumerate
pseudogene biotypes and transcript lengths, retrieve MANE Select IDs, and verify
tool availability. Not part of the pipeline; retained for reproducibility.

**`shared_fraction.py`**  
Reads each `.aln.fa` alignment and slides windows of 75, 100, 150, and 250 bp across
every functional HLA sequence. A window is called "shared" if at least one pseudogene in
the same group has ≥ 90% (or ≥ 95%) nucleotide identity over that window, requiring
at least half the window to be aligned (not gapped) in the pseudogene. Reports the
fraction of windows that are shared per gene, window size, and threshold. Also outputs a
per-position shared flag for plotting (`shared_positions.tsv`). These fractions directly
quantify what proportion of RNA-seq reads from each HLA gene are non-unique by sequence,
as a function of read length — the alignment-space complement of the STAR-based UF
computed in `mappability/sc3UTRs/`.

---

## Shared-sequence fractions (key results)

Computed by `shared_fraction.py` on GENCODE v49 MANE Select vs longest pseudogene
transcripts. Values show % of sliding windows from the functional gene that share
≥ 90% identity with at least one pseudogene in the same group.

| Gene | W=75 bp | W=100 bp | W=150 bp | W=250 bp | Note |
|---|---|---|---|---|---|
| **HLA-A** | 61% | 63% | 62% | 58% | A-clade pseudogenes H, J, L, V |
| **HLA-G** | 58% | 61% | 56% | 52% | |
| HLA-E | 19% | 19% | 15% | 12% | More diverged in this group |
| HLA-B | 22% | 15% | 8% | 4% | BC-clade pseudogenes more diverged |
| HLA-C | 23% | 23% | 11% | 2% | |
| HLA-F | 20% | 16% | 16% | 11% | |
| **HLA-DRB1** | **61%** | **66%** | **72%** | **81%** | DRB pseudogenes 2, 6, 7, 8, 9 |
| **HLA-DRB4** | **62%** | **69%** | **71%** | **79%** | |
| HLA-DRB3 | 56% | 58% | 61% | 70% | |
| HLA-DRB5 | 49% | 50% | 55% | 57% | |
| HLA-DQA1/DQB1 | 0% | 0% | 0% | 0% | DQB3 too diverged |
| HLA-DPA1 | 3% | 2% | 0% | 0% | |
| HLA-DPB1 | 8% | 7% | 6% | 8% | |

**Reading this table:** at a 100 bp read length, 63% of HLA-A reads and 66% of HLA-DRB1
reads share sequence with at least one pseudogene at ≥ 90% identity. These reads cannot
be uniquely attributed to the functional gene. The DRB locus worsens with longer reads
(81% at 250 bp) because the pseudogene coverage, though partial, is concentrated in the
high-conservation exon 2 block that long reads are more likely to fully span. The BC
clade drops sharply with read length, suggesting the pseudogene sequences differ
primarily at short local motifs rather than extended blocks.

At 95% identity threshold (approximate single-nucleotide mismatch tolerance in a 100 bp
window): HLA-A 24%, HLA-DRB1 43%, HLA-DRB4 41% — still far from negligible.

These fractions are **lower bounds** (see caveat above).
