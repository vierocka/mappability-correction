# mappability-correction

In silico mappability correction for RNA-seq quantification at the HLA/APM locus.

This repository documents an ongoing investigation into how tool choice, reference genome configuration, and annotation filtering affect transcript-level quantification in one of the most complex regions of the human genome — the antigen-presenting machinery (APM) locus on chromosome 6. Scripts, figures, and analyses are added incrementally as the project progresses.

---

## Repository structure

```
mappability-correction/
├── mappability/           # Mappability correction scripts (9 parameter variants)
│                          # and supporting documentation — see mappability/README.md
│
├── HLA_MAFFT/             # HLA sequence extraction and MAFFT multiple alignments
│                          # Quantifies sequence identity across HLA paralogs
│
├── Figure_APMlocus/       # Locus map figures (GRCh38.p14) and the script to generate them
│                          # Two-panel annotation overview: TAP/PSM cluster + HLA class I cluster
│
├── Ref/                   # Reference files — NOT tracked (exceed GitHub size limits)
│                          # See ref_files_needed.csv for the file list and download sources
│
├── versions/              # Tool version log
├── ref_files_needed.csv   # Reference file inventory with download sources
└── deduplication_note.txt # CDS deduplication command and stats (seqkit rmdup)
```

---

## Status

Work in progress. Content will be added step by step.

| Module | Status |
|--------|--------|
| APM locus annotation figure | Done |
| HLA sequence identity (MAFFT) | Done |
| Mappability correction scripts | In progress |
| Benchmarking results | Pending |
