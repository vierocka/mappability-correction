#!/usr/bin/env python3
"""
fetch_papers.py — Query OpenAlex for OA human RNA-seq papers in target journals (2020-2026).

Filters applied at query time : journal ISSN + open-access + year range
Filters applied locally        : abstract must contain RNA-seq AND human keywords

Output: papers.tsv
  doi | pmcid | pmid | title | year | journal | n_citations | abstract
"""

import csv, json, re, sys, time
import urllib.request, urllib.parse
from pathlib import Path

HERE  = Path(__file__).resolve().parent
OUT   = HERE / 'papers.tsv'

# OpenAlex polite pool — replace with your email for higher rate limits
EMAIL = 'your@email.com'

# ── Target journals ────────────────────────────────────────────────────────────
# Both print and online ISSNs listed; OpenAlex matches either.
JOURNALS = {
    'Cell Systems':                    ['2405-4712', '2405-4720'],
    'Genome Biology':                  ['1474-760X', '1465-6906'],
    'PLoS Computational Biology':      ['1553-7358', '1553-734X'],
    'Nature Methods':                  ['1548-7105', '1548-7091'],
    'Genome Research':                 ['1088-9051', '1549-5469'],
    'Nature Biotechnology':            ['1087-0156', '1546-1696'],
    'Nature Genetics':                 ['1061-4036', '1546-1718'],
    'eLife':                           ['2050-084X'],
    'Nature Communications':           ['2041-1723'],
    'Nucleic Acids Research':          ['0305-1048', '1362-4962'],
    'Journal of Heredity':             ['0022-1503', '1465-7333'],
    'Bioinformatics':                  ['1367-4803', '1460-2059'],
    'Oncogene':                        ['0950-9232', '1476-5594'],
    'Journal of Clinical Oncology':    ['0732-183X', '1527-7755'],
    'Nature Immunology':               ['1529-2908', '1529-2916'],
    'Scientific Reports':              ['2045-2322'],
    'Molecular Biology and Evolution': ['0737-4038', '1537-1719'],
    'Genome Biology and Evolution':    ['1759-6653'],
}

# ── Local keyword filters ──────────────────────────────────────────────────────
RNA_RE = re.compile(
    r'\b(?:RNA[-\s]?seq|scRNA|snRNA|single[-\s]cell\s+RNA|bulk\s+RNA|'
    r'RNA\s+sequencing|mRNA[-\s]?seq|spatial\s+transcriptom|transcriptom)',
    re.IGNORECASE,
)
HUMAN_RE = re.compile(
    r'\b(?:human|Homo\s+sapiens|patient|donor|cohort|GRCh3[78]|hg(?:19|38)|'
    r'human\s+cell|human\s+tissue|human\s+sample)',
    re.IGNORECASE,
)


def reconstruct_abstract(inv_index):
    """Rebuild plain text from OpenAlex inverted abstract index."""
    if not inv_index:
        return ''
    pos = {}
    for word, positions in inv_index.items():
        for p in positions:
            pos[p] = word
    return ' '.join(pos[i] for i in sorted(pos))


def oa_get(endpoint, params):
    params['mailto'] = EMAIL
    url = endpoint + '?' + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={'User-Agent': f'mailto:{EMAIL}'})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read())


def fetch_journal(name, issns):
    issn_str   = '|'.join(issns)
    base_filter = (f'primary_location.source.issn:{issn_str},'
                   f'is_oa:true,'
                   f'publication_year:2020-2026')
    select = ('id,doi,title,publication_year,cited_by_count,'
              'primary_location,abstract_inverted_index,ids')

    cursor  = '*'
    kept    = []
    page    = 0

    while cursor:
        try:
            data = oa_get('https://api.openalex.org/works', {
                'filter':   base_filter,
                'per-page': 200,
                'cursor':   cursor,
                'select':   select,
            })
        except Exception as exc:
            print(f'  ERROR page {page}: {exc}', file=sys.stderr)
            break

        for w in data.get('results', []):
            title    = w.get('title') or ''
            abstract = reconstruct_abstract(w.get('abstract_inverted_index'))
            text     = title + ' ' + abstract

            if not RNA_RE.search(text):
                continue
            if not HUMAN_RE.search(text):
                continue

            ids   = w.get('ids') or {}
            pmcid_raw = ids.get('pmcid') or ''
            pmcid = re.sub(r'.*/?(PMC\d+)/?$', r'\1', pmcid_raw) if pmcid_raw else ''
            pmid  = re.sub(r'.*/(\d+)/?$', r'\1', ids.get('pmid') or '')
            doi   = (w.get('doi') or '').replace('https://doi.org/', '')

            kept.append({
                'doi':         doi,
                'pmcid':       pmcid,
                'pmid':        pmid,
                'title':       title,
                'year':        w.get('publication_year', ''),
                'journal':     name,
                'n_citations': w.get('cited_by_count', 0),
                'abstract':    abstract,
            })

        cursor = data.get('meta', {}).get('next_cursor')
        page  += 1
        time.sleep(0.12)   # ~8 req/sec — within polite-pool limit

        if page % 10 == 0:
            print(f'    page {page}, {len(kept)} kept')

    return kept


def main():
    all_papers = []
    seen       = set()

    for name, issns in JOURNALS.items():
        print(f'\n{name}')
        papers = fetch_journal(name, issns)
        for p in papers:
            key = p['doi'] or p['pmcid'] or p['title']
            if key and key not in seen:
                seen.add(key)
                all_papers.append(p)
        print(f'  → {len(papers)} matched  (running total: {len(all_papers)})')

    fields = ['doi', 'pmcid', 'pmid', 'title', 'year',
              'journal', 'n_citations', 'abstract']
    with open(OUT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(all_papers)

    print(f'\nSaved {len(all_papers)} papers → {OUT}')


if __name__ == '__main__':
    main()
