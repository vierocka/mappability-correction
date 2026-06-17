#!/usr/bin/env python3
"""
extract_methods.py — Fetch full-text XML from Europe PMC and extract methods sections.

For each paper in papers.tsv:
  1. If PMCID present → Europe PMC fullTextXML → parse <sec>/<app> elements
     whose title matches methods/material/procedure/experimental keywords.
  2. Fallback → Unpaywall PDF URL → pdfplumber text extraction → regex section split.

The script extracts ONLY the methods section, not abstracts or results.
Papers where no methods section is found are recorded with methods_text=''.

Reads:  papers.tsv
Output: methods.tsv  (doi | pmcid | title | year | journal | methods_text | source)
"""

import csv, json, re, sys, time
import urllib.request, urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
IN   = HERE / 'papers.tsv'
OUT  = HERE / 'methods.tsv'

EMAIL = 'your@email.com'   # required for Unpaywall

# When True, skip all network calls and use abstract as the only text source.
# Set to False to attempt Europe PMC + Unpaywall fetches (slow).
ABSTRACT_ONLY = True

# Section titles that indicate a methods section
METHODS_RE = re.compile(
    r'method|material|procedure|experimental\s+design|'
    r'data\s+acqui|subjects?\s+and|patients?\s+and|participants?\s+and|'
    r'clinical\s+data|sample\s+collection|sequencing\s+protocol',
    re.IGNORECASE,
)

MAX_METHODS_CHARS = 12_000   # cap per paper to keep methods.tsv manageable


# ── XML helpers ────────────────────────────────────────────────────────────────

def _tag(el):
    """Return local tag name, stripping XML namespace if present."""
    t = el.tag
    return t.split('}', 1)[-1] if '}' in t else t


def all_text(el):
    """Recursively extract all text from an XML element."""
    parts = [el.text or '']
    for child in el:
        parts.append(all_text(child))
        parts.append(child.tail or '')
    return re.sub(r'\s+', ' ', ' '.join(parts)).strip()


def find_methods_sections(root):
    """
    Walk all <sec> and <app> elements; collect those whose <title> matches
    METHODS_RE. Returns list of text strings (one per matched section).
    """
    hits = []
    for el in root.iter():
        if _tag(el) not in ('sec', 'app'):
            continue
        title_el = el.find('.//{*}title') or el.find('title')
        if title_el is None:
            continue
        title_text = all_text(title_el)
        if METHODS_RE.search(title_text):
            hits.append(all_text(el))
    return hits


# ── Europe PMC ────────────────────────────────────────────────────────────────

def fetch_pmc_xml(pmcid):
    """Return raw XML bytes from Europe PMC, or None on failure."""
    num = re.sub(r'[^0-9]', '', pmcid)   # keep digits only
    url = f'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC{num}/fullTextXML'
    try:
        req = urllib.request.Request(url, headers={'User-Agent': f'mailto:{EMAIL}'})
        with urllib.request.urlopen(req, timeout=40) as r:
            return r.read()
    except Exception:
        return None


def extract_from_pmc(pmcid):
    xml_bytes = fetch_pmc_xml(pmcid)
    if not xml_bytes:
        return None
    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError:
        return None
    sections = find_methods_sections(root)
    if not sections:
        return None
    return '\n\n'.join(sections)[:MAX_METHODS_CHARS]


# ── Unpaywall + PDF fallback ───────────────────────────────────────────────────

def unpaywall_pdf_url(doi):
    url = f'https://api.unpaywall.org/v2/{urllib.parse.quote(doi, safe="")}?email={EMAIL}'
    try:
        with urllib.request.urlopen(url, timeout=20) as r:
            data = json.loads(r.read())
        for loc in (data.get('oa_locations') or []):
            if loc.get('url_for_pdf'):
                return loc['url_for_pdf']
    except Exception:
        pass
    return None


def extract_from_pdf(pdf_url):
    """Download PDF and extract methods text via pdfplumber + regex."""
    try:
        import pdfplumber, io
    except ImportError:
        return None
    try:
        req = urllib.request.Request(
            pdf_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=60) as r:
            raw = r.read()
        pages = []
        with pdfplumber.open(io.BytesIO(raw)) as pdf:
            for page in pdf.pages:
                pages.append(page.extract_text() or '')
        text = '\n'.join(pages)
        # Regex-based section extraction from plain text
        m = re.search(
            r'(?:^|\n)\s*'
            r'(?:Materials?\s+and\s+Methods?|Methods?\s+and\s+Materials?|'
            r'Methods?|Experimental\s+Procedures?|STAR\s*Methods?|'
            r'Patients?\s+and\s+Methods?|Subjects?\s+and\s+Methods?)\s*\n'
            r'(.+?)'
            r'(?=\n\s*(?:Results?|Discussion|Conclusion|Acknowledgements?|'
            r'References?|Supplementary)\s*\n|\Z)',
            text, re.IGNORECASE | re.DOTALL,
        )
        if m:
            return m.group(1)[:MAX_METHODS_CHARS]
    except Exception:
        pass
    return None


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    papers = []
    with open(IN) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            papers.append(row)
    print(f'Loaded {len(papers)} papers')

    # Build doi → abstract lookup for fallback
    abstracts = {p.get('doi', ''): (p.get('abstract') or '').strip()
                 for p in papers}

    fields = ['doi', 'pmcid', 'title', 'year', 'journal', 'methods_text', 'source']
    n_pmc = n_pdf = n_abstract = n_fail = 0

    with open(OUT, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t',
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader()

        for i, paper in enumerate(papers):
            pmcid = (paper.get('pmcid') or '').strip()
            doi   = (paper.get('doi')   or '').strip()
            methods_text = ''
            source       = 'not_found'

            if not ABSTRACT_ONLY:
                # ── Path 1: Europe PMC full-text XML ──────────────────────────
                if pmcid:
                    methods_text = extract_from_pmc(pmcid) or ''
                    if methods_text:
                        source = 'europepmc_xml'
                        n_pmc += 1
                    time.sleep(0.15)

                # ── Path 2: Unpaywall → PDF ────────────────────────────────────
                if not methods_text and doi:
                    pdf_url = unpaywall_pdf_url(doi)
                    if pdf_url:
                        methods_text = extract_from_pdf(pdf_url) or ''
                        if methods_text:
                            source = 'unpaywall_pdf'
                            n_pdf += 1
                    time.sleep(0.3)

            # ── Path 3: Abstract fallback (always available) ───────────────────
            if not methods_text:
                ab = abstracts.get(doi, '')
                if ab:
                    methods_text = ab
                    source = 'abstract'
                    n_abstract += 1
                else:
                    n_fail += 1

            writer.writerow({
                'doi':          doi,
                'pmcid':        pmcid,
                'title':        paper.get('title', ''),
                'year':         paper.get('year', ''),
                'journal':      paper.get('journal', ''),
                'methods_text': methods_text,
                'source':       source,
            })

            if (i + 1) % 500 == 0:
                pct = 100 * (i + 1) / len(papers)
                print(f'  {i+1}/{len(papers)} ({pct:.0f}%)  '
                      f'PMC={n_pmc}  PDF={n_pdf}  abstract={n_abstract}  failed={n_fail}',
                      flush=True)

    print(f'\nDone.')
    print(f'  Europe PMC XML : {n_pmc}')
    print(f'  PDF fallback   : {n_pdf}')
    print(f'  Abstract       : {n_abstract}')
    print(f'  Not found      : {n_fail}')
    print(f'  Saved → {OUT}')


if __name__ == '__main__':
    main()
