# Comparison scheme — mappability correction project

Overview of what is compared, why, and the escalation path when short-read RNA-seq is insufficient.

---

## Comparison design

| Read source | STAR reference | Kallisto index | Expected bias | Script series |
|---|---|---|---|---|
| MANE RNA | Full genome (705 seqs) | — | ⬇ Signal LOSS — multimapper discard | v1, v3, v4, v6 |
| MANE RNA | Dedup CDS (93,088 seqs) | — | ⬇ Partial LOSS — fewer competitors | v2 |
| MANE RNA | MANE RNA (19,437 seqs) | — | UF → 1.0 — self-mapping upper bound | v5 |
| MANE RNA | — | MANE index (UTRs) | ⬆ INFLATION — pseudogene collapse + UTR rescue | Kallisto crossval |
| MANE RNA | — | Dedup CDS index | ⬆ INFLATION moderate — no UTR rescue | Kallisto crossval |
| Exon-flank | Full genome / dedup CDS | — | Intronic reads — intron retention proxy | v1–v2 (flank) |
| Exon-flank | — | Both indices | Unmapped fraction — expected pseudoalignment failure | Kallisto crossval |
| Dedup CDS | Full genome | — | CDS genomic uniqueness — identifies problematic coding loci | v7_a–f, v7_i–j |
| Dedup CDS | CDS self-map | — | Within-CDS redundancy — residual overlap after dedup | v7_g–h |
| Dedup CDS | — | Both indices | CDS representation in pseudoalignment space | Kallisto crossval |

UTR smoothing contribution = (Kallisto TPM with MANE index) − (Kallisto TPM with dedup CDS index), for the same MANE RNA read set.

---

## Flowchart — bias, classification, and escalation

```mermaid
flowchart TD

    classDef loss fill:#fca5a5,stroke:#dc2626,color:#111
    classDef inflation fill:#fed7aa,stroke:#ea580c,color:#111
    classDef neutral fill:#e0e7ff,stroke:#4338ca,color:#111
    classDef ok fill:#bbf7d0,stroke:#16a34a,color:#111
    classDef warn fill:#fef9c3,stroke:#b45309,color:#111
    classDef fail fill:#fecaca,stroke:#b91c1c,color:#111
    classDef metric fill:#bfdbfe,stroke:#1d4ed8,color:#111
    classDef best fill:#ddd6fe,stroke:#6d28d9,color:#111

    %%──── READ SOURCES ──────────────────────────────────────────────────────────
    R1["🧬 MANE spliced RNA reads
    19,437 transcripts · SE+PE · 75–200 bp
    polyA library proxy"]
    R2["🔬 Exon-flank reads
    exon ±50 bp intron · SE+PE · 75–150 bp
    intron retention proxy"]
    R3["📦 Dedup CDS reads  (v7)
    93,088 CDS sequences · SE+PE · 75–200 bp
    coding-region uniqueness"]

    %%──── STAR ALIGNMENTS ───────────────────────────────────────────────────────
    subgraph STAR_BOX ["STAR + featureCounts  ➜  uniqueness factor  UF = MAPQ-255 fraction"]
        direction LR
        ST_G["v1 / v3 / v4 / v6
        STAR → full genome
        705 seqs"]:::loss
        ST_C["v2
        STAR → dedup CDS
        93,088 seqs"]:::loss
        ST_M["v5
        STAR → MANE RNA
        19,437 seqs"]:::neutral
        ST_V["v7_a–f · v7_i–j
        STAR → full genome
        CDS reads"]:::metric
        ST_S["v7_g–h
        STAR → CDS self-map
        CDS reads"]:::metric
    end

    R1 --> ST_G & ST_C & ST_M
    R2 --> ST_G & ST_C
    R3 --> ST_V & ST_S

    %%──── KALLISTO ──────────────────────────────────────────────────────────────
    subgraph KALL_BOX ["Kallisto  ➜  TPM per transcript  ➜  inflation estimate"]
        direction LR
        KM["MANE RNA reads
        → MANE index
        UTRs included"]:::inflation
        KC["MANE RNA reads
        → dedup CDS index
        no UTRs"]:::inflation
        KF["Exon-flank reads
        → both indices
        expected: fails pseudoalignment"]:::fail
        KV["CDS reads
        → both indices
        CDS representation in pseudoalignment"]:::metric
    end

    R1 --> KM & KC
    R2 --> KF
    R3 --> KV

    %%──── BIAS OUTCOMES ─────────────────────────────────────────────────────────
    ST_G -->|"multimapper discard
    typical UF 0.4 – 0.9
    gene-body loss worst"| B_LOSS["⬇ Signal LOSS"]:::loss

    ST_C -->|"fewer competing loci
    UF higher than full genome
    but not zero"| B_LOSS

    ST_M -->|"no pseudogene competition
    UF → 1.0 upper bound"| B_MANE["MANE self-map baseline
    UF ≈ 1.0"]:::neutral

    ST_V -->|"CDS reads vs genomic decoys
    UF per CDS sequence
    identifies unquantifiable coding loci"| B_CDS["CDS genomic uniqueness
    independent of transcript annotation"]:::metric

    ST_S -->|"residual k-mer overlap
    after exact-duplicate removal"| B_CDS

    KM -->|"pseudogene reads collapse +
    UTR k-mer rescue
    highest inflation"| B_INF["⬆ Signal INFLATION"]:::inflation

    KC -->|"pseudogene collapse
    no UTR rescue
    moderate inflation"| B_INF

    KF -->|"intronic reads
    fail pseudoalignment
    unmapped signal"| B_IR["Unmapped fraction
    ≈ intron retention
    in real tumour data"]:::fail

    KV --> B_CDS

    %%──── METRICS ───────────────────────────────────────────────────────────────
    B_LOSS & B_MANE --> UF[["Uniqueness factor
    UF = MAPQ-255 reads / simulated reads
    per MANE transcript
    reference-specific"]]:::metric

    B_INF --> INF[["Inflation estimate
    ratio TPM_MANE_idx / TPM_CDS_idx
    via Kallisto crossval
    isolates UTR smoothing component"]]:::metric

    B_CDS --> UF_CDS[["CDS uniqueness factor  (v7)
    UF per CDS sequence
    identifies coding loci with poor
    genomic uniqueness"]]:::metric

    %%──── GENE CLASSIFICATION ───────────────────────────────────────────────────
    UF & INF & UF_CDS --> CLS{"Gene classification
    per MANE transcript"}

    CLS -->|"UF ≈ 1.0
    no inflation"| G1["✅  Always OK
    any pipeline reliable
    single-copy genes, unique loci"]:::ok

    CLS -->|"0 < UF < 1
    no or mild inflation"| G2["🔧  Apply correction
    true ≈ observed / UF
    RPL/RPS partial, KZNF, protocadherins
    many cancer-testis antigens"]:::warn

    CLS -->|"UF ≈ 1.0
    inflated vs CDS index"| G3["⚠  Prefer STAR + fC full genome
    Kallisto overestimates
    pseudogenes absent from MANE
    not transcribed in this context"]:::warn

    CLS -->|"UF ≈ 0
    and inflated"| G4["❌  Always wrong
    both tools biased in opposite directions
    cannot cancel
    ─────────────────────────────
    HLA-A/B/C with full genome ref
    canonical histones  ·  ubiquitin
    RPL/RPS high-density pseudogenes
    olfactory receptors
    HERV-K expressed loci"]:::fail

    %%──── ESCALATION ────────────────────────────────────────────────────────────
    G4 --> PROT(["Escalation
    Protein-level measurement
    proteomics / immunopeptidomics"])

    PROT -->|"unique tryptic peptides
    correct allele in database
    no pseudogene translation conflict"| POK["✅  Protein level reliable
    allele-specific quantification
    when database is complete"]:::ok

    PROT --> PF1["❌  Shared peptides between paralogs
    HLA-A/B/C: conserved α3 domain
    RPL/RPS: most peptides non-unique
    ubiquitin: identical CDS
    ──────────────────────────────
    protein inference ambiguity
    = RNA multimapping in peptide space
    Nesvizhskii & Aebersold  MCP 2005"]:::fail

    PROT --> PF2["❌  Allele absent from search database
    HLA α1/α2 peptides unmatched
    → allele invisible in results
    engine choice shifts IDs by 30–70%
    ──────────────────────────────
    Parker et al.  MCP 2021
    Pak et al.  MCP 2021"]:::fail

    PROT --> PF3["❌  Database too large / FDR distortion
    all isoforms + variants + non-canonical
    → FDR anti-conservative
    true hits suppressed by decoy inflation
    ──────────────────────────────
    Li et al.  BMC Genomics 2016"]:::fail

    PROT --> PF4["❌  Pseudogene translation & non-canonical ORFs
    ≥1,546 pseudogenes produce detectable peptides
    non-canonical proteins populate immunopeptidome
    post-translational spliced peptides in HLA
    ──────────────────────────────
    Vasylieva  JPR 2024
    Ruiz Cuevas et al.  Cell Rep 2021
    Mishto et al.  Proteomics 2022"]:::fail

    %%──── BEST SOLUTION ─────────────────────────────────────────────────────────
    PF1 & PF2 & PF3 & PF4 --> LONG["Long reads + HLA typing +
    allele-matched immunopeptidomics
    ─────────────────────────────────────
    PacBio / ONT + WGS HLA typing
    personalised protein reference DB
    DIA-MS + sample-matched spectral library
    ─────────────────────────────────────
    Pak et al. MCP 2021 — multi-HLA spectral libraries
    Marcu et al. JITC 2021 — HLA Ligand Atlas"]:::best

    LONG --> BEST(["Best achievable resolution
    allele-specific quantification
    ─────────────────────────────────
    Still limited by:
    · pseudogene-derived peptides not in any DB
    · post-translationally spliced peptides
    · somatic HLA mutations absent from reference
    · HLA loss of heterozygosity in tumour"]):::best
```

---

## Summary

The central result is that **no single tool and reference combination is unbiased for all genes**. Signal loss (STAR + full genome) and signal inflation (Kallisto + MANE) are not random noise — they are systematic, locus-specific, and directional. Any prediction model, immune score, or cross-cohort comparison that uses RNA-seq expression from pseudogene-dense loci without correction inherits this bias as a fixed systematic error.

Protein-level measurement resolves some cases but faces structurally identical problems in peptide space: shared peptides (= RNA multimapping), missing alleles in the database (= reference mismatch), and pseudogene-derived peptides absent from any database (= Kallisto inflation from non-canonical sources).

The most accurate allele-specific quantification requires long reads with HLA typing and a sample-matched immunopeptidomics database — and even then, somatic HLA mutations and post-translational peptide splicing remain unresolved in standard workflows.
