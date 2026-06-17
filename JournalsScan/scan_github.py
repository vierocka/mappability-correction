#!/usr/bin/env python3
"""
scan_github.py — Fetch pipeline scripts from GitHub repos found in tools_summary.tsv
                 and extract alignment/quantification parameters.

Reads:   tools_summary.tsv  (github_links column from mine_tools.py)

Output:  github_scan.tsv   — one row per repo with extracted pipeline parameters

Authentication:
  Set GITHUB_TOKEN env variable for 5000 req/hour.
  Without it: 60 req/hour (fine for small repo sets, slow for large ones).

  export GITHUB_TOKEN=ghp_...
"""

import base64, csv, json, os, re, sys, time
import urllib.request, urllib.error
from pathlib import Path
from collections import defaultdict

HERE   = Path(__file__).resolve().parent
IN     = HERE / 'tools_summary.tsv'
OUT    = HERE / 'github_scan.tsv'

GITHUB_TOKEN = os.environ.get('GITHUB_TOKEN', '')
API          = 'https://api.github.com'

# Files worth fetching from a repo
WORKFLOW_RE = re.compile(
    r'(?:Snakefile|.*\.smk|.*\.sh|.*\.nf|nextflow\.config'
    r'|.*pipeline.*\.(?:py|sh|r|R|yaml|yml)'
    r'|.*workflow.*\.(?:py|sh|r|R)'
    r'|.*align.*\.(?:sh|py)|.*mapping.*\.(?:sh|py)'
    r'|.*run.*\.sh|config\.(?:yaml|yml|json)'
    r'|README(?:\.md|\.rst|\.txt)?)',
    re.IGNORECASE,
)

# Parameters to extract from pipeline files
EXTRACT = {
    'star_genome_fa':       re.compile(r'--genomeFastaFiles\s+(\S+)'),
    'star_gtf':             re.compile(r'--sjdbGTFfile\s+(\S+)'),
    'star_genome_dir':      re.compile(r'--genomeDir\s+(\S+)'),
    'star_multNmax':        re.compile(r'--outSAMmultNmax\s+(-?\d+)'),
    'star_multimapNmax':    re.compile(r'--outFilterMultimapNmax\s+(\d+)'),
    'star_soft_clip_none':  re.compile(r'--alignEndsType\s+EndToEnd'),
    'featurecounts_gtf':    re.compile(r'featureCounts\s+[^#\n]*-a\s+(\S+)'),
    'featurecounts_M':      re.compile(r'featureCounts\s+[^#\n]*\s-M\b'),
    'featurecounts_O':      re.compile(r'featureCounts\s+[^#\n]*\s-O\b'),
    'featurecounts_s':      re.compile(r'featureCounts\s+[^#\n]*-s\s+([012])'),
    'salmon_decoys':        re.compile(r'--decoys\s+(\S+)'),
    'salmon_index_fasta':   re.compile(r'salmon\s+index\s+[^#\n]*(?:-t|--transcripts)\s+(\S+)'),
    'kallisto_index_fasta': re.compile(r'kallisto\s+index\s+[^#\n]*\s(\S+\.fa(?:\.gz)?)'),
    'cellranger_fasta':     re.compile(r'cellranger\s+mkref[^#\n]*--fasta\s+(\S+)'),
    'cellranger_genes':     re.compile(r'cellranger\s+mkref[^#\n]*--genes\s+(\S+)'),
    'gencode_version':      re.compile(r'gencode[._\-]?v?(\d+)', re.IGNORECASE),
    'ensembl_version':      re.compile(r'(?:ensembl|release)[._\s](\d{2,3})', re.IGNORECASE),
    'grch38':               re.compile(r'GRCh38(?:\.p\d+)?|hg38'),
    'grch37':               re.compile(r'GRCh37|hg19'),
    'primary_assembly':     re.compile(r'primary[_\s]assembly|no[_\s]alt', re.IGNORECASE),
}


def github_request(path, params=None):
    url = API + path
    if params:
        url += '?' + '&'.join(f'{k}={v}' for k, v in params.items())
    headers = {'User-Agent': 'mappability-scan/1.0', 'Accept': 'application/vnd.github+json'}
    if GITHUB_TOKEN:
        headers['Authorization'] = f'Bearer {GITHUB_TOKEN}'
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            remaining = int(r.headers.get('X-RateLimit-Remaining', 999))
            if remaining < 5:
                reset_ts = int(r.headers.get('X-RateLimit-Reset', 0))
                wait = max(0, reset_ts - time.time()) + 2
                print(f'  Rate limit low ({remaining} remaining). Waiting {wait:.0f}s …')
                time.sleep(wait)
            return json.loads(r.read())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None   # repo not found / private
        if e.code == 403:
            print(f'  403 Forbidden for {path} — rate limited or token issue')
            time.sleep(60)
            return None
        raise
    except Exception:
        return None


def get_default_branch(owner, repo):
    data = github_request(f'/repos/{owner}/{repo}')
    if data is None:
        return None
    return data.get('default_branch', 'main')


def list_repo_files(owner, repo, branch):
    """Return list of (path, url) for workflow-relevant files using the Trees API."""
    data = github_request(
        f'/repos/{owner}/{repo}/git/trees/{branch}',
        {'recursive': '1'}
    )
    if data is None or 'tree' not in data:
        return []
    files = []
    for item in data['tree']:
        if item.get('type') != 'blob':
            continue
        p = item['path']
        if WORKFLOW_RE.search(p.split('/')[-1]):
            files.append((p, item.get('url', '')))
    return files


def fetch_file_content(blob_url):
    """Fetch file content from GitHub blob URL (base64 encoded)."""
    if not blob_url:
        return None
    data = github_request(blob_url.replace(API, ''))
    if data is None:
        return None
    content_b64 = data.get('content', '')
    if not content_b64:
        return None
    try:
        return base64.b64decode(content_b64).decode('utf-8', errors='replace')
    except Exception:
        return None


def scan_text(text):
    """Apply all EXTRACT patterns to text, return dict of first matches."""
    results = {}
    for key, pat in EXTRACT.items():
        m = pat.search(text)
        if m:
            results[key] = m.group(1) if m.lastindex else 'True'
    return results


def scan_repo(owner_repo):
    """Scan one GitHub repo. Returns dict of extracted parameters."""
    parts = owner_repo.replace('github.com/', '').split('/')
    if len(parts) < 2:
        return {'error': 'bad_url'}
    owner, repo = parts[0], parts[1]

    branch = get_default_branch(owner, repo)
    if branch is None:
        return {'error': 'not_found_or_private'}
    time.sleep(0.3)

    files = list_repo_files(owner, repo, branch)
    time.sleep(0.3)

    if not files:
        return {'error': 'no_workflow_files', 'branch': branch}

    combined = {}
    files_scanned = []
    for fpath, blob_url in files[:30]:   # cap at 30 files per repo
        content = fetch_file_content(blob_url)
        time.sleep(0.2)
        if not content:
            continue
        hits = scan_text(content)
        if hits:
            files_scanned.append(fpath)
            for k, v in hits.items():
                if k not in combined:   # first match wins
                    combined[k] = v

    combined['files_scanned'] = ';'.join(files_scanned)
    combined['n_files_scanned'] = len(files_scanned)
    combined['branch'] = branch
    return combined


def infer_multimapper_policy(scan):
    """
    Classify how multimappers are handled based on extracted parameters.
    Returns: 'likely_discard' | 'likely_keep_all' | 'em_redistribute' | 'explicit_keep' | 'unknown'
    """
    multNmax = scan.get('star_multNmax')
    if multNmax is not None:
        if multNmax == '-1':
            return 'explicit_keep_all'
        if multNmax == '1':
            return 'unique_only'
        if int(multNmax) > 1:
            return 'keep_top_N'

    # featureCounts -M means count multimappers
    if scan.get('featurecounts_M'):
        return 'explicit_keep'

    # EM-based tools always redistribute
    if scan.get('salmon_index_fasta') or scan.get('salmon_decoys'):
        return 'em_redistribute'
    if scan.get('kallisto_index_fasta'):
        return 'em_redistribute'

    # STAR + featureCounts without -M → multimappers discarded by featureCounts
    if scan.get('star_genome_dir') and scan.get('featurecounts_gtf'):
        return 'likely_discard'

    return 'unknown'


def main():
    # Collect unique github links from tools_summary.tsv
    repo_to_papers = defaultdict(list)
    with open(IN) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            links = row.get('github_links', '')
            if not links:
                continue
            doi = row.get('doi', '')
            for link in links.split(';'):
                link = link.strip()
                if link:
                    repo_to_papers[link].append(doi)

    repos = sorted(repo_to_papers)
    print(f'Found {len(repos)} unique GitHub repos across {sum(len(v) for v in repo_to_papers.values())} paper links')
    if not GITHUB_TOKEN:
        print('WARNING: GITHUB_TOKEN not set. Rate limit = 60 req/hour. '
              'Set export GITHUB_TOKEN=... for faster scanning.')

    out_fields = [
        'repo', 'doi_list', 'error', 'branch', 'n_files_scanned', 'files_scanned',
        'multimapper_policy',
        'star_genome_fa', 'star_gtf', 'star_genome_dir',
        'star_multNmax', 'star_multimapNmax', 'star_soft_clip_none',
        'featurecounts_gtf', 'featurecounts_M', 'featurecounts_O', 'featurecounts_s',
        'salmon_decoys', 'salmon_index_fasta', 'kallisto_index_fasta',
        'cellranger_fasta', 'cellranger_genes',
        'gencode_version', 'ensembl_version', 'grch38', 'grch37', 'primary_assembly',
    ]

    with open(OUT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields, delimiter='\t',
                           extrasaction='ignore', restval='')
        w.writeheader()

        for i, repo in enumerate(repos):
            print(f'  [{i+1}/{len(repos)}] {repo}', end=' ', flush=True)
            scan = scan_repo(repo)
            policy = infer_multimapper_policy(scan)
            print(f'→ {policy}  ({scan.get("n_files_scanned", 0)} files)')

            row = {
                'repo':               repo,
                'doi_list':           ';'.join(repo_to_papers[repo]),
                'error':              scan.get('error', ''),
                'multimapper_policy': policy,
            }
            row.update(scan)
            w.writerow(row)

    print(f'\nSaved → {OUT}')


if __name__ == '__main__':
    main()
