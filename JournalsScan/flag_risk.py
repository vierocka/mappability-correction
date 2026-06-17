#!/usr/bin/env python3
"""
flag_risk.py — Flag papers at risk of mappability-driven misinterpretation.

The core concern: studies of HLA/MHC genes or immune-related loci that use
DE signatures to classify patients may be affected by systematic mappability
artifacts. The problem is *differential*: cancer-specific pseudogene
transcription changes the multimapping landscape between tumor and normal,
so the bias is not constant and cannot cancel out in DE analysis.

Risk logic
----------
CRITICAL : HLA/MHC in methods + patient classification + tool discards multimappers
           → apparent HLA downregulation may be misread as immune evasion

HIGH     : HLA/MHC in methods + patient classification + EM tool without decoys
           → reads redistributed across paralogs/pseudogenes → inflated or
             deflated allele-level signal; DE signatures unreliable

MEDIUM   : HLA/MHC OR pseudogene mention + DE signature + cancer context
           → risk depends on tool choice and reference completeness

LOW      : Any of the above but tool handling verified (explicit -M, decoys used,
           multNmax=-1) OR no patient classification

Reads:
  tools_summary.tsv   (from mine_tools.py)
  github_scan.tsv     (from scan_github.py) — optional; improves accuracy

Output:
  risk_flags.tsv      — one row per paper, risk level + rationale
  risk_summary.tsv    — counts by risk level, journal, year
"""

import csv
from collections import Counter, defaultdict
from pathlib import Path

HERE      = Path(__file__).resolve().parent
TOOLS_IN  = HERE / 'tools_summary.tsv'
GITHUB_IN = HERE / 'github_scan.tsv'
OUT_RISK  = HERE / 'risk_flags.tsv'
OUT_SUM   = HERE / 'risk_summary.tsv'


def truthy(val):
    return str(val).strip().lower() in ('true', '1', 'yes')


def load_github_policies():
    """Load multimapper policy per repo from github_scan.tsv (if it exists)."""
    if not GITHUB_IN.exists():
        return {}
    policies = {}
    with open(GITHUB_IN) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            repo = row.get('repo', '')
            if repo:
                policies[repo] = row.get('multimapper_policy', 'unknown')
    return policies


def github_policy_for_paper(github_links, repo_policies):
    """Return the verified multimapper policy if any repo for this paper was scanned."""
    for link in (github_links or '').split(';'):
        link = link.strip()
        if link in repo_policies:
            return repo_policies[link]
    return None


def assess_risk(row, verified_policy):
    """
    Returns (risk_level, rationale_list).
    risk_level: 'CRITICAL' | 'HIGH' | 'MEDIUM' | 'LOW' | 'NONE'
    """
    hla          = truthy(row.get('bio_HLA'))
    mhc          = truthy(row.get('bio_MHC'))
    pseudo       = truthy(row.get('bio_pseudogene'))
    cancer       = truthy(row.get('bio_cancer'))
    de_sig       = truthy(row.get('bio_DE_signature'))
    classify_pt  = truthy(row.get('bio_classify_patient'))
    immune_ev    = truthy(row.get('bio_immune_evasion'))
    immunotherapy= truthy(row.get('bio_immunotherapy'))
    amp          = truthy(row.get('bio_AMP'))
    multi_note   = truthy(row.get('bio_multimapper_note'))

    uses_star    = truthy(row.get('aln_STAR'))
    uses_hisat   = truthy(row.get('aln_HISAT2'))
    uses_cr      = truthy(row.get('aln_CellRanger'))
    uses_starsolo= truthy(row.get('aln_STARsolo'))
    uses_fcounts = truthy(row.get('cnt_featureCounts'))
    uses_htseq   = truthy(row.get('cnt_HTSeq'))
    uses_kallisto= truthy(row.get('aln_kallisto'))
    uses_salmon  = truthy(row.get('aln_Salmon'))
    uses_rsem    = truthy(row.get('cnt_RSEM'))
    uses_alevin  = truthy(row.get('aln_Alevin'))
    uses_kb      = truthy(row.get('aln_kbpython'))

    hla_context  = hla or mhc
    em_tool      = uses_kallisto or uses_salmon or uses_rsem or uses_alevin or uses_kb
    discard_tool = ((uses_star or uses_hisat or uses_cr or uses_starsolo) and
                    (uses_fcounts or uses_htseq))
    any_aligner  = (uses_star or uses_hisat or uses_cr or uses_starsolo or
                    uses_kallisto or uses_salmon or uses_rsem or uses_alevin or uses_kb)

    rationale = []

    # Check if GitHub scan confirmed explicit multimapper handling
    explicitly_handled = False
    if verified_policy in ('explicit_keep_all', 'explicit_keep', 'em_redistribute'):
        explicitly_handled = True

    # ── CRITICAL ──────────────────────────────────────────────────────────────
    if hla_context and classify_pt and discard_tool and not explicitly_handled:
        rationale.append(
            'HLA/MHC studied + patient classified + tool discards multimappers '
            '(STAR+featureCounts or HISAT2+HTSeq without -M); apparent HLA expression '
            'changes may reflect pseudogene-driven multimapper loss, not biology.'
        )
        if immune_ev or immunotherapy:
            rationale.append(
                'Paper also discusses immune evasion or immunotherapy: '
                'misattributing pipeline artifact as HLA loss is clinically dangerous.'
            )
        return 'CRITICAL', rationale

    # ── HIGH ──────────────────────────────────────────────────────────────────
    if hla_context and classify_pt and em_tool and not explicitly_handled:
        rationale.append(
            'HLA/MHC studied + patient classified + EM tool (kallisto/Salmon/RSEM/Alevin); '
            'reads ambiguous between HLA alleles and pseudogenes are redistributed by EM, '
            'inflating or deflating allele-specific signal. If pseudogene reference not included, '
            'signal accumulates on canonical alleles → biased DE signature.'
        )
        return 'HIGH', rationale

    if amp and classify_pt and any_aligner and not explicitly_handled:
        rationale.append(
            'AMP cohort/locus mentioned + patient classified; '
            'AMP RA/SLE studies frequently use HLA-associated DE signatures — '
            'mappability handling should be verified.'
        )
        return 'HIGH', rationale

    # ── MEDIUM ────────────────────────────────────────────────────────────────
    if (hla_context or pseudo) and de_sig and cancer:
        rationale.append(
            'HLA/MHC or pseudogene mentioned + DE signature used + cancer context; '
            'cancer-specific pseudogene derepression changes multimapping landscape '
            'between tumor and normal → differential bias not cancelled by DE.'
        )
        if explicitly_handled:
            rationale.append(f'GitHub scan: multimapper policy = {verified_policy}.')
        return 'MEDIUM', rationale

    if hla_context and de_sig and not classify_pt:
        rationale.append(
            'HLA/MHC in DE analysis; patient-level classification not detected, '
            'but HLA DE results may still be affected by multimapper handling.'
        )
        return 'MEDIUM', rationale

    # ── LOW ───────────────────────────────────────────────────────────────────
    if (hla_context or pseudo or amp) and any_aligner:
        rationale.append(
            'HLA/MHC, pseudogene, or AMP context detected; no patient classification '
            'or DE signature extracted from methods. Risk is present but limited.'
        )
        if multi_note:
            rationale.append('Paper mentions multimapping — may be aware of the issue.')
        return 'LOW', rationale

    return 'NONE', []


def main():
    repo_policies = load_github_policies()
    n_repos = len(repo_policies)
    print(f'Loaded {n_repos} repo policies from github_scan.tsv'
          if n_repos else 'github_scan.tsv not found — using methods text only')

    rows = []
    with open(TOOLS_IN) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            rows.append(row)
    print(f'Loaded {len(rows)} papers')

    risk_counts = Counter()
    year_risk   = defaultdict(Counter)
    journal_risk= defaultdict(Counter)
    out_rows    = []

    for row in rows:
        github_links   = row.get('github_links', '')
        verified_policy = github_policy_for_paper(github_links, repo_policies)
        risk, rationale = assess_risk(row, verified_policy)

        risk_counts[risk] += 1
        year    = row.get('year', '')
        journal = row.get('journal', '')
        if year:    year_risk[year][risk] += 1
        if journal: journal_risk[journal][risk] += 1

        out_rows.append({
            'doi':              row.get('doi', ''),
            'pmcid':            row.get('pmcid', ''),
            'title':            row.get('title', ''),
            'year':             year,
            'journal':          journal,
            'n_citations':      row.get('n_citations', ''),
            'risk_level':       risk,
            'rationale':        ' | '.join(rationale),
            'verified_policy':  verified_policy or '',
            'combo':            row.get('combo', ''),   # present if combo_summary ran first
            'bio_HLA':          row.get('bio_HLA', ''),
            'bio_MHC':          row.get('bio_MHC', ''),
            'bio_pseudogene':   row.get('bio_pseudogene', ''),
            'bio_cancer':       row.get('bio_cancer', ''),
            'bio_DE_signature': row.get('bio_DE_signature', ''),
            'bio_classify_patient': row.get('bio_classify_patient', ''),
            'bio_immune_evasion':   row.get('bio_immune_evasion', ''),
            'bio_immunotherapy':    row.get('bio_immunotherapy', ''),
            'bio_AMP':          row.get('bio_AMP', ''),
            'bio_multimapper_note': row.get('bio_multimapper_note', ''),
            'likely_scrna':     row.get('likely_scrna', ''),
            'likely_bulk':      row.get('likely_bulk', ''),
            'github_links':     github_links,
        })

    # Sort: CRITICAL first, then HIGH, MEDIUM, LOW, NONE
    RANK = {'CRITICAL': 0, 'HIGH': 1, 'MEDIUM': 2, 'LOW': 3, 'NONE': 4}
    out_rows.sort(key=lambda r: (RANK.get(r['risk_level'], 5),
                                  -(int(r['n_citations']) if r['n_citations'] else 0)))

    # ── Write risk_flags.tsv ─────────────────────────────────────────────────
    fields = [
        'risk_level', 'doi', 'pmcid', 'title', 'year', 'journal', 'n_citations',
        'rationale', 'verified_policy', 'combo',
        'bio_HLA', 'bio_MHC', 'bio_pseudogene', 'bio_cancer',
        'bio_DE_signature', 'bio_classify_patient',
        'bio_immune_evasion', 'bio_immunotherapy', 'bio_AMP', 'bio_multimapper_note',
        'likely_scrna', 'likely_bulk', 'github_links',
    ]
    with open(OUT_RISK, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t',
                           extrasaction='ignore', restval='')
        w.writeheader()
        w.writerows(out_rows)
    print(f'Saved {len(out_rows)} rows → {OUT_RISK}')

    # ── Write risk_summary.tsv ────────────────────────────────────────────────
    levels  = ['CRITICAL', 'HIGH', 'MEDIUM', 'LOW', 'NONE']
    years   = sorted(year_risk)
    journals= sorted(journal_risk)
    n_total = len(rows)

    sum_rows = []
    for lvl in levels:
        r = {
            'risk_level': lvl,
            'n_papers':   risk_counts[lvl],
            'pct':        f'{100*risk_counts[lvl]/n_total:.1f}' if n_total else '0',
        }
        for yr in years:
            r[yr] = year_risk[yr].get(lvl, 0)
        sum_rows.append(r)

    sum_fields = ['risk_level', 'n_papers', 'pct'] + years
    with open(OUT_SUM, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=sum_fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(sum_rows)
    print(f'Saved → {OUT_SUM}')

    # ── Console summary ───────────────────────────────────────────────────────
    print(f'\nRisk distribution (n={n_total} papers):')
    print(f'\n{"Level":<10} {"n":>6} {"pct":>6}   ' + '  '.join(f'{y}' for y in years))
    print('─' * (30 + 6 * len(years)))
    for lvl in levels:
        n   = risk_counts[lvl]
        pct = f'{100*n/n_total:.1f}%' if n_total else '0%'
        yr_cols = '  '.join(f'{year_risk[y].get(lvl,0):>4}' for y in years)
        print(f'{lvl:<10} {n:>6} {pct:>6}   {yr_cols}')

    # Per-journal breakdown for CRITICAL + HIGH
    print('\nCRITICAL + HIGH papers by journal:')
    for jrn in sorted(journal_risk, key=lambda j: -(journal_risk[j]['CRITICAL'] +
                                                     journal_risk[j]['HIGH'])):
        nc = journal_risk[jrn]['CRITICAL']
        nh = journal_risk[jrn]['HIGH']
        if nc + nh == 0:
            continue
        print(f'  {jrn:<45} CRITICAL={nc}  HIGH={nh}')


if __name__ == '__main__':
    main()
