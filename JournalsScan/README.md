# JournalsScan

Automated pipeline to identify published RNA-seq studies where **mappability artifacts may have shaped gene expression results** — specifically studies that use multi-mapping-sensitive tools (e.g. kallisto/Salmon EM redistribution) or that discard multi-mappers without correction and then draw conclusions in genomic regions known to be affected (HLA locus, pseudogene-rich regions, antisense overlaps).

## Motivation

The mappability correction analysis (parent folder) shows that ~7.5% of MANE transcripts have non-unique 3' ends under 10x Genomics-style 100 bp PE sequencing, and that HLA genes, tandem gene families, and genes with antisense overlaps are disproportionately affected.  A natural question is: **which published studies could have inflated or deflated expression calls for these genes, and is the effect large enough to change patient classification, DE signatures, or cell-type annotation?**

Rather than re-analysing raw data, this pipeline screens the methods sections of papers to:

1. Identify the aligner + counter combination (which determines whether multi-mappers are discarded, redistributed, or collapsed)
2. Flag papers that work in affected biological domains (HLA typing, immune evasion, cancer classification, scRNA-seq cell annotation)
3. Cross-reference GitHub repositories for analysis code that could be inspected or re-run

The output is not a list of retracted papers but a **prioritised reading list** for follow-up: papers where the combination of tool choice and biological domain makes mappability-driven signal inflation most plausible.

## Pipeline

```
fetch_papers.py          Download paper metadata + abstracts from PubMed / Europe PMC
        ↓
extract_methods.py       Extract methods text (full-text XML if OA, else abstract)
        ↓
upgrade_methods_fulltext.py   Upgrade abstract rows to full-text for OA papers (incremental)
        ↓
mine_tools.py            Detect tools, reference genomes, annotation versions, GitHub links
        ↓
scan_github.py           Inspect linked GitHub repos for analysis scripts / data files
        ↓
combo_summary.py         Summarise tool combinations and assign risk categories
        ↓
flag_risk.py             Per-paper risk flags + summary table
```

## Risk categories

| Risk | Meaning |
|------|---------|
| `em_redistribute` | kallisto or Salmon used — multi-mappers redistributed by EM; expression may be inflated for poorly unique genes |
| `discard_default` | STAR / CellRanger / featureCounts with default settings — multi-mappers discarded; expression deflated for low-UF genes |
| `unknown` | Tool combination detected but no clear mappability consequence, or no tool detected |

Risk is escalated to **CRITICAL/HIGH** when a tool with a known mappability effect is combined with HLA / immune-evasion biology or patient classification.

## Key outputs

| File | Contents |
|------|---------|
| `papers.tsv` | Source paper metadata (DOI, PMID, title, year, journal, abstract) |
| `methods.tsv` | Methods text per paper (source: `europepmc_xml`, `abstract`, or `not_found`) |
| `tools_summary.tsv` | Per-paper tool flags, annotation versions, GitHub links, snippets |
| `aggregate_summary.tsv` | Tool prevalence counts broken down by year |
| `combo_summary.tsv` | Unique tool combinations with risk label and study-type breakdown |
| `risk_flags.tsv` | Per-paper risk level with reasons |
| `risk_summary.tsv` | Distribution of risk levels by year |

## Configuration

`extract_methods.py` has an `ABSTRACT_ONLY = True` flag at the top.  Set it to `False` to attempt Europe PMC full-text and Unpaywall PDF fetching (slower, requires network access).

`upgrade_methods_fulltext.py` upgrades rows that have a PMID but only an abstract; it resolves PMCIDs via Europe PMC and fetches full-text XML for Open Access papers.  Run it after `extract_methods.py` to improve tool detection from ~1.5% (abstract-only) toward the expected ~30–50% (full-text).

Set `EMAIL` at the top of each script to a valid address for polite API access (Unpaywall and Europe PMC both request it).

## Limitations

- Abstract-only mode detects biological terms well but tool names poorly (~1.5% of papers show a tool hit vs ~30–50% expected from full-text).
- GitHub scanning is rate-limited to 60 requests/hour without a `GITHUB_TOKEN`; set `export GITHUB_TOKEN=...` before running `scan_github.py`.
- The risk categories are heuristic; a paper flagged `em_redistribute` may still report robust results if the affected genes are not in the conclusions.
