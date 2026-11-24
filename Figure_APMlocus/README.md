# Figure_APMlocus — annotation map of the HLA/APM locus

Two-panel locus overview of the antigen processing machinery (APM) and HLA class I region on chromosome 6p21.3 (GRCh38.p14). The figure makes visible the annotation geometry that underlies the multimapping problem in RNA-seq.

## Why this locus is important

### The antigen processing machinery

The APM is the cellular system that generates and loads peptides onto HLA class I molecules for presentation at the cell surface. Key components in this chromosomal region include:

- **HLA-A, HLA-B, HLA-C** — the classical HLA class I genes encoding the heavy chain of the peptide-MHC complex.
- **PSMB8 (LMP7) and PSMB9 (LMP2)** — immunoproteasome subunits that replace constitutive subunits upon IFN-γ stimulation and shift the cleavage specificity toward peptides with hydrophobic C-termini, which bind HLA class I groove more efficiently.
- **TAP1 and TAP2** — the transporter associated with antigen processing, which translocates proteasome-generated peptides from the cytosol into the ER lumen.
- **TAPBP (tapasin)** — the scaffold of the peptide loading complex (PLC), bridging TAP and HLA class I / beta-2-microglobulin (B2M) and facilitating high-affinity peptide selection.

B2M (the invariant beta-2-microglobulin subunit) is encoded on chromosome 15 and is not shown in this figure, but it is required for stable HLA class I surface expression.

### Complex annotation geometry

The HLA/APM region on 6p21.3 is one of the most gene-dense and polymorphic regions in the human genome. What makes it structurally unusual from an RNA-seq standpoint:

1. **Pseudogenes.** Multiple HLA pseudogenes (e.g., HLA-H, HLA-J, HLA-K, HLA-L) and APM-related pseudogenes are embedded within the region, sharing high sequence identity with functional paralogs.
2. **Alternative loci.** GRCh38 includes numerous alt loci for HLA genes (over 300 alt sequences across the MHC region), representing common haplotypic variants. Reads from a sample carrying a non-reference allele may align to both the primary assembly and one or more alt loci simultaneously.
3. **Overlapping gene bodies.** Some genes in this region have overlapping transcripts or nested annotation, producing reads that cannot be unambiguously assigned to a single gene even with perfect mapping.
4. **Extreme polymorphism.** HLA-A, -B, and -C are the most polymorphic human genes (thousands of alleles each at the protein level). Reads from alleles divergent from GRCh38 may fail to map uniquely.

### Relevance to immunopeptidomics

Immunopeptidomics — the mass spectrometry-based identification of peptides presented on HLA class I or class II molecules — depends directly on the correct expression and function of all APM components. Quantification errors in any of these genes propagate into:

- Inaccurate modelling of peptide presentation probability (used by tools such as NetMHCpan and MHCflurry).
- Misclassification of tumours as antigen-presentation competent or deficient.
- Incorrect RNA-to-protein correlation at the locus when benchmarking immunopeptidomics against transcriptomics.

The same allele-mismatch and shared-peptide problems that limit RNA-seq at this locus also affect mass spectrometry-based HLA peptide identification — the two technologies share a structural ambiguity problem in different spaces (k-mer vs peptide-spectrum).

### Relevance to immuno-oncology

HLA class I expression on tumour cells is a prerequisite for cytotoxic T-cell recognition. Its deregulation is a primary immune evasion mechanism:

- **HLA loss of heterozygosity (LOH)** occurs in 15–40% of many solid tumour types and causes complete or partial loss of HLA allele surface expression, preventing T-cell recognition of tumour antigens.
- **APM downregulation** — silencing of PSMB8, PSMB9, TAP1, TAP2, or TAPBP impairs the peptide loading pathway without HLA gene loss, producing a functional antigen presentation deficiency invisible to HLA copy-number analysis.
- **IFN-γ pathway mutations** (JAK1, JAK2, B2M) abolish the HLA upregulation response to immune signalling, conferring resistance to checkpoint immunotherapy.

Immunotherapy response prediction tools and antigen presentation scores (TIS, IPS, TIDE, CYT) depend on accurate RNA-seq quantification of exactly these genes. Multimapping bias at this locus directly corrupts model inputs, and its magnitude varies by tool and reference in a systematic, correctable way — which is what this project characterises.

## Figure files

| File | Description |
|---|---|
| `APM_locus_chr6.pdf` | Main figure (vector, publication-quality) |
| `APM_locus_chr6.png` | Raster version for quick viewing |
| `APM_locus_GFF_snapshot.tsv` | Annotation snapshot from MANE GFF used to build the figure |
| `plot_APM_locus.py` | Script to reproduce the figure from the GFF snapshot |
