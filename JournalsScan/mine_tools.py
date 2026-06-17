#!/usr/bin/env python3
"""
mine_tools.py — Extract aligners, reference genomes, annotation versions,
                counters, and GitHub links from methods sections.

Reads:  methods.tsv  (from extract_methods.py)
        papers.tsv   (for citation counts)

Output:
  tools_summary.tsv    — one row per paper, all extracted fields
  aggregate_summary.tsv — counts per tool across all papers, broken down by year
"""

import csv, re, sys
from collections import defaultdict, Counter
from pathlib import Path

HERE          = Path(__file__).resolve().parent
METHODS_TSV   = HERE / 'methods.tsv'
PAPERS_TSV    = HERE / 'papers.tsv'
TOOLS_OUT     = HERE / 'tools_summary.tsv'
AGGREGATE_OUT = HERE / 'aggregate_summary.tsv'

CTX = 130   # characters of surrounding context saved per matched snippet


# ── Patterns ──────────────────────────────────────────────────────────────────

PATTERNS = {

    # ── Reference genome version ──────────────────────────────────────────────
    'ref_GRCh38':   re.compile(r'\b(?:GRCh38(?:\.p\d+)?|hg38)\b'),
    'ref_GRCh37':   re.compile(r'\b(?:GRCh37|hg19)\b'),
    'ref_T2T':      re.compile(r'\bT2T[-_]CHM13\b|Telomere.to.Telomere', re.IGNORECASE),

    # ── Annotation source ─────────────────────────────────────────────────────
    'ann_GENCODE':  re.compile(r'\bGENCODE\b', re.IGNORECASE),
    'ann_Ensembl':  re.compile(r'\bEnsembl\b',  re.IGNORECASE),
    'ann_RefSeq':   re.compile(r'\bRefSeq\b|GCF_\d{9}\.\d+'),
    'ann_UCSC':     re.compile(r'\bUCSC\b(?:\s+genome|\s+browser|\s+annotation)',
                               re.IGNORECASE),

    # ── Reference processing ──────────────────────────────────────────────────
    'proc_ERCC':      re.compile(r'\bERCC\b'),
    'proc_spike':     re.compile(r'\bspike.in\b', re.IGNORECASE),
    'proc_primary':   re.compile(r'primary\s+(?:assembly|chromosome|contigs?)',
                                 re.IGNORECASE),
    'proc_MANE':      re.compile(r'\bMANE\b'),
    'proc_sjdb':      re.compile(r'\bsjdbOverhang\b|\bsjdb\b', re.IGNORECASE),
    'proc_filter_mt': re.compile(r'(?:remov|filter|exclud).{0,40}mitochondri',
                                 re.IGNORECASE),

    # ── Bulk RNA aligners ─────────────────────────────────────────────────────
    # Negative lookahead avoids matching "STAR Methods" section header
    'aln_STAR':     re.compile(r'\bSTAR\b(?!\s*[✦\*]?\s*[Mm]ethod)'),
    'aln_HISAT2':   re.compile(r'\bHISAT2\b',   re.IGNORECASE),
    'aln_TopHat':   re.compile(r'\bTopHat[12]?\b', re.IGNORECASE),
    'aln_bowtie':   re.compile(r'\bbowtie[12]?\b', re.IGNORECASE),
    'aln_subread':  re.compile(r'\bsubread\b',   re.IGNORECASE),

    # ── Pseudoaligners / quasi-mappers ────────────────────────────────────────
    'aln_kallisto': re.compile(r'\bkallisto\b',  re.IGNORECASE),
    # exclude "Salmon" as a fish or food item
    'aln_Salmon':   re.compile(
        r'\bSalmon\b(?!\s+(?:protein|fish|trout|spawn|dish|fillet))', re.IGNORECASE),

    # ── scRNA aligners / preprocessors ───────────────────────────────────────
    'aln_CellRanger': re.compile(r'\bCell\s?Ranger\b',  re.IGNORECASE),
    'aln_STARsolo':   re.compile(r'\bSTARsolo\b',       re.IGNORECASE),
    'aln_Alevin':     re.compile(r'\bAlevin(?:-fry)?\b', re.IGNORECASE),
    'aln_kbpython':   re.compile(
        r'\bkb[-_]python\b|\bkb\s+ref\b|\bkallisto[|/\\]bustools\b', re.IGNORECASE),
    'aln_zUMIs':      re.compile(r'\bzUMIs\b'),
    'aln_Optimus':    re.compile(r'\bOptimus\b(?:\s+pipeline)?'),  # HCA pipeline

    # ── Counters / quantifiers ────────────────────────────────────────────────
    'cnt_featureCounts': re.compile(r'\bfeatureCounts\b', re.IGNORECASE),
    'cnt_HTSeq':         re.compile(r'\bHTSeq\b'),
    'cnt_RSEM':          re.compile(r'\bRSEM\b'),
    'cnt_bustools':      re.compile(r'\bbustools\b',      re.IGNORECASE),
    'cnt_simpleaf':      re.compile(r'\bsimpleaf\b|\bAlevin-fry\b', re.IGNORECASE),
    'cnt_StringTie':     re.compile(r'\bStringTie\b',     re.IGNORECASE),

    # ── Downstream tools (confirm RNA-seq workflow type) ──────────────────────
    'ds_DESeq2':   re.compile(r'\bDESeq2\b'),
    'ds_edgeR':    re.compile(r'\bedgeR\b'),
    'ds_limma':    re.compile(r'\blimma\b',    re.IGNORECASE),
    'ds_Seurat':   re.compile(r'\bSeurat\b'),
    'ds_Scanpy':   re.compile(r'\bScanpy\b',   re.IGNORECASE),
    'ds_Monocle':  re.compile(r'\bMonocle[23]?\b'),
    'ds_Harmony':  re.compile(r'\bHarmony\b'),
    'ds_scVI':     re.compile(r'\bscVI\b|\bSCVI\b'),

    # ── Biological domain terms ───────────────────────────────────────────────
    # These flag studies where mappability artifacts in HLA / pseudogene-rich
    # regions may have shaped DE signatures used for patient classification.
    'bio_HLA':              re.compile(
        r'\bHLA[-\s]?[A-Z][A-Z0-9]*\b|human\s+leukocyte\s+antigen',
        re.IGNORECASE),
    'bio_MHC':              re.compile(
        r'\bMHC(?:\s+class\s+[I]{1,3})?\b|major\s+histocompatibility',
        re.IGNORECASE),
    'bio_pseudogene':       re.compile(r'\bpseudogene\b', re.IGNORECASE),
    'bio_AMP':              re.compile(
        r'\bAMP\s+(?:locus|RA|SLE|lupus|cohort|project|phase|network)\b'
        r'|Accelerated\s+Med(?:icine)?s?\s+Partnership',
        re.IGNORECASE),
    'bio_DE_signature':     re.compile(
        r'(?:transcriptional|gene|expression)\s+signature'
        r'|differentially\s+expressed\s+gene|DEGs?\b'
        r'|gene\s+set\s+(?:enrichment|score)',
        re.IGNORECASE),
    'bio_classify_patient': re.compile(
        r'(?:classif\w+|stratif\w+|predict\w+)\s+.{0,80}'
        r'(?:patient|sample|disease\s+subtype|clinical\s+outcome|responder|survival)',
        re.IGNORECASE),
    'bio_immune_evasion':   re.compile(
        r'immune\s+(?:evasion|escape|surveillance)'
        r'|MHC.{0,40}(?:loss|downregul|reduc|absent)'
        r'|HLA.{0,40}(?:loss|downregul|absent|silenc)',
        re.IGNORECASE),
    'bio_immunotherapy':    re.compile(
        r'\bimmunotherap\w+\b|checkpoint\s+inhibit\w+'
        r'|anti[-\s]?(?:PD[-\s]?[L1]+|CTLA[-\s]?4)\b|CAR[-\s]?T\b',
        re.IGNORECASE),
    'bio_cancer':           re.compile(
        r'\b(?:cancer|tumor|tumour|malignant|carcinoma|melanoma|lymphoma'
        r'|leukemia|sarcoma|glioma|adenocarcinoma|neoplasm)\b',
        re.IGNORECASE),
    'bio_multimapper_note': re.compile(
        r'multi.?mapp\w+|multiple.{0,30}mapping|NH:i:'
        r'|outSAMmultNmax|outFilterMultimapNmax',
        re.IGNORECASE),
}

GENCODE_VER_RE = re.compile(
    r'\bGENCODE\s+(?:release\s+|version\s+)?v?(\d+)\b', re.IGNORECASE)
ENSEMBL_VER_RE = re.compile(
    r'\bEnsembl\s+(?:release\s+|version\s+)?(\d+)\b', re.IGNORECASE)
GITHUB_RE = re.compile(
    r'github\.com/([\w\-][\w\-\.]*?/[\w\-][\w\-\.]*?)(?=[,\s;:)\]>\'\"<]|$)',
    re.IGNORECASE)


# ── Helpers ───────────────────────────────────────────────────────────────────

def snippet(text, m):
    s = max(0, m.start() - CTX)
    e = min(len(text), m.end() + CTX)
    return '...' + text[s:e].replace('\n', ' ').replace('\t', ' ') + '...'


def clean_github(owner_repo):
    """Remove trailing punctuation and .git from a captured owner/repo string."""
    u = re.sub(r'[.,;:)\]>\'\"]+$', '', owner_repo)
    u = re.sub(r'\.git$', '', u, flags=re.IGNORECASE)
    return 'github.com/' + u


def mine(text):
    hits     = {}
    snippets = {}

    for key, pat in PATTERNS.items():
        m = pat.search(text)
        hits[key]     = bool(m)
        snippets[key] = snippet(text, m) if m else ''

    # Versioned annotations
    hits['gencode_versions'] = ';'.join(sorted(set(GENCODE_VER_RE.findall(text))))
    hits['ensembl_versions'] = ';'.join(sorted(set(ENSEMBL_VER_RE.findall(text))))

    # GitHub links — deduplicate, clean trailing punctuation
    raw_links = GITHUB_RE.findall(text)
    hits['github_links'] = ';'.join(sorted(set(clean_github(u) for u in raw_links)))

    # Inferred study type
    hits['likely_scrna'] = any(hits.get(k) for k in [
        'aln_CellRanger', 'aln_STARsolo', 'aln_Alevin', 'aln_kbpython',
        'aln_zUMIs', 'aln_Optimus', 'ds_Seurat', 'ds_Scanpy', 'ds_Monocle',
        'ds_Harmony', 'ds_scVI',
    ])
    hits['likely_bulk'] = any(hits.get(k) for k in [
        'aln_STAR', 'aln_HISAT2', 'aln_TopHat', 'aln_Salmon', 'aln_kallisto',
        'cnt_featureCounts', 'cnt_HTSeq', 'cnt_RSEM', 'cnt_StringTie',
        'ds_DESeq2', 'ds_edgeR', 'ds_limma',
    ])

    return hits, snippets


# ── Output field list ─────────────────────────────────────────────────────────

BOOL_FIELDS = list(PATTERNS.keys()) + ['likely_scrna', 'likely_bulk']
STR_FIELDS  = ['gencode_versions', 'ensembl_versions', 'github_links']
META_FIELDS = ['doi', 'pmcid', 'title', 'year', 'journal', 'n_citations', 'source']
SNIP_FIELDS = [f'snip_{k}' for k in PATTERNS]
ALL_FIELDS  = META_FIELDS + BOOL_FIELDS + STR_FIELDS + SNIP_FIELDS


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    # Load citation counts from papers.tsv
    cites = {}
    with open(PAPERS_TSV) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            cites[row['doi']] = row.get('n_citations', '')

    # Process methods
    rows_out = []
    tool_counts  = Counter()   # total count per tool
    year_counts  = defaultdict(Counter)   # year → tool → count

    n_ok = n_fail = 0

    with open(METHODS_TSV) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            mtext = (row.get('methods_text') or '').strip()

            if not mtext:
                n_fail += 1
                continue
            n_ok += 1

            hits, snippets = mine(mtext)
            year = row.get('year', '')
            doi  = row.get('doi', '')

            out = {f: row.get(f, '') for f in META_FIELDS}
            out['n_citations'] = cites.get(doi, '')
            out.update({k: hits[k] for k in BOOL_FIELDS + STR_FIELDS})
            out.update({f'snip_{k}': snippets.get(k, '') for k in PATTERNS})
            rows_out.append(out)

            for k in BOOL_FIELDS:
                if hits.get(k):
                    tool_counts[k] += 1
                    if year:
                        year_counts[year][k] += 1

    # ── Write tools_summary.tsv ───────────────────────────────────────────────
    with open(TOOLS_OUT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=ALL_FIELDS, delimiter='\t',
                           extrasaction='ignore', restval='')
        w.writeheader()
        w.writerows(rows_out)
    print(f'Saved {len(rows_out)} rows → {TOOLS_OUT}')

    # ── Write aggregate_summary.tsv ───────────────────────────────────────────
    years = sorted(year_counts)
    agg_fields = ['tool', 'n_papers', 'pct'] + years
    n_total = len(rows_out)

    with open(AGGREGATE_OUT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=agg_fields, delimiter='\t')
        w.writeheader()
        for tool, n in sorted(tool_counts.items(), key=lambda x: -x[1]):
            row = {
                'tool':     tool,
                'n_papers': n,
                'pct':      f'{100*n/n_total:.1f}' if n_total else '0',
            }
            for yr in years:
                row[yr] = year_counts[yr].get(tool, 0)
            w.writerow(row)
    print(f'Saved → {AGGREGATE_OUT}')

    # ── Console summary ───────────────────────────────────────────────────────
    print(f'\nPapers mined: {n_ok}  |  no methods found: {n_fail}')
    print(f'\n{"Tool":<25} {"n":>6} {"pct":>6}   ' + '  '.join(f'{y}' for y in years))
    print('─' * (40 + 6 * len(years)))
    for tool, n in sorted(tool_counts.items(), key=lambda x: -x[1]):
        pct = f'{100*n/n_total:.1f}%' if n_total else '0%'
        yr_cols = '  '.join(f'{year_counts[y].get(tool, 0):>4}' for y in years)
        print(f'{tool:<25} {n:>6} {pct:>6}   {yr_cols}')


if __name__ == '__main__':
    main()
