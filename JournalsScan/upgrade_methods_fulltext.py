#!/usr/bin/env python3
"""
upgrade_methods_fulltext.py  —  Upgrade methods.tsv from abstract to full-text
                                 for Open Access papers via Europe PMC.

For each row in methods.tsv with source='abstract' that has a PMID:
  1. Query Europe PMC search API to resolve PMCID from PMID.
  2. If PMCID found, fetch fullTextXML and extract methods section.
  3. Replace abstract with methods text, update source to 'europepmc_xml'.

Papers without PMID or not in PMC stay as abstract.

Reads:  methods.tsv, papers.tsv (for PMID lookup)
Output: methods.tsv  (in-place upgrade, overwrites)

Usage:
    python3 upgrade_methods_fulltext.py
"""

import csv, json, re, time, sys
import urllib.request, urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
METHODS = HERE / 'methods.tsv'
PAPERS  = HERE / 'papers.tsv'

EMAIL       = 'your@email.com'
SLEEP_SRCH  = 0.18   # seconds between Europe PMC search calls
SLEEP_XML   = 0.18   # seconds between full-text XML calls
MAX_CHARS   = 12_000

METHODS_RE = re.compile(
    r'method|material|procedure|experimental\s+design|'
    r'data\s+acqui|subjects?\s+and|patients?\s+and|participants?\s+and|'
    r'clinical\s+data|sample\s+collection|sequencing\s+protocol',
    re.IGNORECASE,
)


def _tag(el):
    t = el.tag
    return t.split('}', 1)[-1] if '}' in t else t


def all_text(el):
    parts = [el.text or '']
    for child in el:
        parts.append(all_text(child))
        parts.append(child.tail or '')
    return re.sub(r'\s+', ' ', ' '.join(parts)).strip()


def find_methods_sections(root):
    hits = []
    for el in root.iter():
        if _tag(el) not in ('sec', 'app'):
            continue
        title_el = el.find('.//{*}title') or el.find('title')
        if title_el is None:
            continue
        if METHODS_RE.search(all_text(title_el)):
            hits.append(all_text(el))
    return hits


def pmid_to_pmcid(pmid):
    """Return PMCID (e.g. 'PMC12345') or None."""
    url = (f'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
           f'?query=EXT_ID:{pmid}%20AND%20SRC:MED'
           f'&resultType=lite&format=json')
    try:
        req = urllib.request.Request(url, headers={'User-Agent': f'mailto:{EMAIL}'})
        with urllib.request.urlopen(req, timeout=20) as r:
            data = json.loads(r.read())
        for result in data.get('resultList', {}).get('result', []):
            pmcid = result.get('pmcid', '')
            if pmcid:
                return pmcid
    except Exception:
        pass
    return None


def fetch_fulltext(pmcid):
    """Return extracted methods text or None."""
    num = re.sub(r'[^0-9]', '', pmcid)
    url = f'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC{num}/fullTextXML'
    try:
        req = urllib.request.Request(url, headers={'User-Agent': f'mailto:{EMAIL}'})
        with urllib.request.urlopen(req, timeout=40) as r:
            xml_bytes = r.read()
        root = ET.fromstring(xml_bytes)
        sections = find_methods_sections(root)
        if not sections:
            return None
        return '\n\n'.join(sections)[:MAX_CHARS]
    except Exception:
        return None


def main():
    # Load PMID lookup: doi → pmid
    doi_to_pmid = {}
    with open(PAPERS) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            pmid = (row.get('pmid') or '').strip()
            doi  = (row.get('doi') or '').strip()
            if pmid and doi:
                doi_to_pmid[doi] = pmid

    # Load current methods.tsv
    rows = []
    with open(METHODS) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    print(f'Loaded {len(rows)} rows from methods.tsv')

    fields = list(rows[0].keys()) if rows else \
             ['doi', 'pmcid', 'title', 'year', 'journal', 'methods_text', 'source']

    candidates = [r for r in rows
                  if r.get('source') == 'abstract'
                  and doi_to_pmid.get(r.get('doi', ''))]
    print(f'Candidates for upgrade (abstract + has PMID): {len(candidates)}')

    upgraded = 0
    pmcid_found = 0
    n_done = 0

    # Build a dict for fast update
    row_by_doi = {r['doi']: r for r in rows}

    for i, row in enumerate(candidates):
        doi  = row.get('doi', '')
        pmid = doi_to_pmid.get(doi, '')

        # Step 1: resolve PMCID
        pmcid = pmid_to_pmcid(pmid)
        time.sleep(SLEEP_SRCH)

        if not pmcid:
            n_done += 1
            if (i + 1) % 500 == 0:
                print(f'  {i+1}/{len(candidates)} searched  '
                      f'pmcid_found={pmcid_found}  upgraded={upgraded}', flush=True)
            continue

        pmcid_found += 1
        row_by_doi[doi]['pmcid'] = pmcid

        # Step 2: fetch full-text
        methods_text = fetch_fulltext(pmcid)
        time.sleep(SLEEP_XML)

        if methods_text:
            row_by_doi[doi]['methods_text'] = methods_text
            row_by_doi[doi]['source']       = 'europepmc_xml'
            upgraded += 1

        n_done += 1
        if n_done % 200 == 0:
            print(f'  {n_done}/{len(candidates)} searched  '
                  f'pmcid_found={pmcid_found}  upgraded={upgraded}', flush=True)

    # Write updated methods.tsv
    print(f'\nUpgraded {upgraded} rows from abstract → europepmc_xml')
    print(f'PMCID resolved: {pmcid_found}/{len(candidates)}')
    with open(METHODS, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t',
                           lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    print(f'Saved → {METHODS}')


if __name__ == '__main__':
    main()
