#!/usr/bin/env python3
"""
Apply the agreed uniqueness / multimapping rules to raw featureCounts output,
then average corrected TPM across replicates within each condition.

Inputs per sample, in <cell>/raw_counts/ (2026-10-07):
  <SRR>_<cell>_MANEpersUTR_fcM.counts.csv  featureCounts -M --fraction   -> COUNTS
  <SRR>_<cell>_MANEpersUTR_fc.counts.csv   featureCounts, unique only    -> fc for the
                                            fc ~ fcM check and HLA allele ratios
  <SRR>_<cell>_mane_fc.counts.csv          basic MANE, unique only       -> HLA-A/B/C totals
MANEpersUTR = basic MANE with the sample's own HLA-A/B/C alleles (incl. UTRs); every
other transcript ID and length is identical to basic MANE.

Rules (decided 2026-10-07):
  1. Counts = fcM for every transcript (fractional multimappers included).
  2. Alt-haplotype duplicates whose primary twin is also in MANE (FOLD_INTO_PRIMARY,
     same gene / same protein): the alt fcM is added to the primary copy, the alt row
     is excluded.
  3. Uniqueness factor f applied ONLY to pseudogene inflation: genome factor < 0.9,
     pseudogenes = TRUE in low_uniqueness_MANE_transcripts.tsv, and no visible
     competition inside MANE in that sample (shared = 1 - fc/fcM < SHARED_MAX, or
     fcM < MIN_FCM_FOR_SHARED). f = 0 -> ZERO_FLOOR.
  4. f = 1 for everything else: factor >= 0.9, paralogs / duplicates in MANE with real
     sharing (fcM > fc), self-homology, nothing detected, re-included non-primary.
  5. HLA-A/B/C: total = basic-MANE canonical count; f by rules 3-4 like every other
     gene (HLA-A: UF < 0.9 + pseudogenes -> x UF; HLA-B/-C: f = 1). Heterozygous
     lines (HCC1937, HCT116): split into the two allele rows (IDs as in MANEpersUTR,
     e.g. ENST00000376809.10_A0201) by the ratio of UNIQUE allele reads in
     MANEpersUTR_fc of that sample. JY (homozygous) keeps the canonical row.
  6. Excluded (TPM_corrected = NaN, out of the normalisation): non-primary
     excluded=TRUE (except CELL_SPECIFIC_INCLUDE), the 17 PAR genes, non-coding
     transcripts (non_coding_MANE_transcripts_no_protein.tsv).

  counts_corrected = (counts [+ folded alt]) x f_unique          (never divided)
  TPM_corrected    = RPK_corrected / sum(RPK_corrected) x 1e6    (within the sample)

Factor table per cell line = real library read length (STAR logs): HCC1937 L75,
HCT116 / JY L100.

Conditions:
  HCC1937: 74_75   -- SRR31477874 + SRR31477875 (SRR925709 = older batch, not mixed)
  HCT116:  41_42   -- SRR3644141 + SRR3644142
  JY:      364065  -- SRR364065

Outputs, per cell line, in <cell>/corrected_counts/:
  <SRR>_mane_fc_corrected.tsv                              per sample (+ `rule`, `shared`)
  <cell>_corrected_mean_<condition>_mane_fc_corrected.tsv  per condition mean

Usage:
    python3 apply_uniqueness_and_average.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent

# --- PIPELINE HOOK (2026-11-22): base dir + config from Snakemake ---------
import json as _json, os as _os
BASE = Path(_os.environ.get("PIPE_BASE", BASE))
_CFG = _json.load(open(_os.environ["PIPE_CONFIG"])) if _os.environ.get("PIPE_CONFIG") else {}
# -------------------------------------------------------------------------

FACTORS_DIR = BASE / "uniqueness_factors"
LOW_UF_FILE = FACTORS_DIR / "low_uniqueness_MANE_transcripts.tsv"
EXCLUDED_FILE = BASE / "non_primary_MANE_transcripts.tsv"   # 64 non-primary MANE transcripts; `excluded` column decides
NONCODING_FILE = BASE / "non_coding_MANE_transcripts_no_protein.tsv"

ZERO_FLOOR = 0.001        # f == 0 -> retain 0.1% of reads (only reached via rule 3)
UF_MAX = 0.9              # rule 3 only below this genome factor; f = 1 at or above
SHARED_MAX = 0.10         # rule 3: shared = 1 - fc/fcM below this = no MANE competition
MIN_FCM_FOR_SHARED = 20   # below this many fcM reads the shared fraction is not used

HLA = {"ENST00000376809.10": "HLA-A", "ENST00000412585.7": "HLA-B",
       "ENST00000376228.10": "HLA-C"}

# Rule 2: alt-haplotype copy -> primary-chromosome copy of the same gene, both in MANE
FOLD_INTO_PRIMARY = {
    "ENST00000421274.7": "ENST00000428956.7",   # C4A (alt)       -> C4A
    "ENST00000435500.7": "ENST00000435363.7",   # C4B_2 (alt)     -> C4B
    "ENST00000615140.2": "ENST00000619989.1",   # CCL3L3 (alt)    -> CCL3L3
    "ENST00000634014.2": "ENST00000620576.5",   # CCL4L2 (alt)    -> CCL4L2
    "ENST00000621406.5": "ENST00000380729.8",   # GTF2H2C_2 (alt) -> GTF2H2C
    "ENST00000570508.3": "ENST00000249116.7",   # APOBEC3A (alt)  -> APOBEC3A
    "ENST00000618279.3": "ENST00000290765.9",   # GSTT2 (alt)     -> GSTT2B
}

# Pseudoautosomal genes: identical X/Y copies (UF = 0), not quantifiable per copy
PAR_GENES = {
    "ENST00000313871.9": "AKAP17A", "ENST00000381241.9": "ASMT", "ENST00000381317.9": "ASMTL",
    "ENST00000381192.10": "CD99", "ENST00000400841.8": "CRLF2", "ENST00000381529.9": "CSF2RA",
    "ENST00000334651.11": "DHRSX", "ENST00000326153.10": "GTPBP6", "ENST00000331035.10": "IL3RA",
    "ENST00000244174.11": "IL9R", "ENST00000381297.10": "P2RY8", "ENST00000381657.8": "PLCXD1",
    "ENST00000390665.9": "PPP2R3B", "ENST00000686671.1": "SHOX", "ENST00000381401.11": "SLC25A6",
    "ENST00000286448.12": "VAMP7", "ENST00000652001.1": "ZBED1",
}

# Fix-patch MANE transcripts brought back in (2026-10-07). Their MANE locus is a
# GRC fix patch (= corrected reference sequence, not a population variant), so the
# primary-assembly simulation never produced a factor (table value empty, status
# excluded_fix_patch). Checked: each maps to exactly ONE primary-assembly locus
# (minimap2 -x splice, MAPQ 60, 93-100 % identity = the patch correction) with no
# competing copy at >=30 % coverage / >=90 % identity except minor pseudogene hits
# (IFNL4P1, OR4C4P, a processed ZNF598 pseudogene on chr9) that leave the primary
# locus at MAPQ 60; the same genes are single-copy primary genes in Ensembl 116
# (nf-core standard). Computed uniqueness = unique placement -> f_unique = 1.0.
# Still excluded from the fix-patch group: NBPF1 (no unique locus, 7 hits at
# MAPQ 0-1), PRAMEF22 and ENSG00000292227 (absent from the primary annotation).
INCLUDED_FIX_PATCH_F = 1.0
INCLUDED_FIX_PATCH = {
    "ENST00000644422.3": "ABO",
    "ENST00000710274.1": "ASPN",
    "ENST00000671703.3": "CASP8AP2",
    "ENST00000709217.2": "CD300H",
    "ENST00000628705.3": "FCGBP",
    "ENST00000643531.4": "FOXO6",
    "ENST00000709220.2": "GBA3",
    "ENST00000710324.2": "GPATCH4",
    "ENST00000710321.1": "IFNL4",
    "ENST00000707626.2": "MROH5",
    "ENST00000710289.2": "MROH8",
    "ENST00000711671.1": "MUC16",
    "ENST00000713550.1": "MUC2",
    "ENST00000710306.2": "NKRF",
    "ENST00000710315.1": "OR10J4",
    "ENST00000708727.2": "OR11H7",
    "ENST00000710319.1": "OR4C45",
    "ENST00000710368.2": "OR4Q2",
    "ENST00000710310.1": "OR52E1",
    "ENST00000710330.1": "OR5AL1",
    "ENST00000715685.1": "ORAI1",
    "ENST00000709157.1": "PNLIPRP2",
    "ENST00000643490.2": "POLR2A",
    "ENST00000640199.4": "RFLNB",
    "ENST00000643872.4": "RYBP",
    "ENST00000673075.2": "SAMD1",
    "ENST00000710353.1": "SHANK3",
    "ENST00000642844.3": "SLC37A4",
    "ENST00000708883.3": "TCEAL6",
    "ENST00000713554.1": "VN1R5",
    "ENST00000710288.1": "ZNF598",
}

# Alt-haplotype transcripts whose alt copy is MANE's ONLY representative (the
# primary twin is not in MANE), MAPQ 60, no extra genomic loci (2026-10-07).
INCLUDED_ALT_SINGLE = {
    "ENST00000709635.1": "FUT3",
    "ENST00000709446.2": "CASP12",
    "ENST00000703726.1": "OR12D1",
}

# All transcripts with excluded=FALSE in non_primary_MANE_transcripts.tsv, f = 1.0
INCLUDED_ALL = {**INCLUDED_FIX_PATCH, **INCLUDED_ALT_SINGLE}

# Cell-line-specific inclusions of otherwise excluded alt-haplotype genes, f = 1.0.
# HLA-DRB3 / HLA-DRB4 exist only on DR52 / DR53 haplotypes (column
# DR52_DR53_haplotype); JY carries and expresses both (13,988 / 11,394 raw reads,
# MAPQ 60, 92-94 % identity to HLA-DRB1).
CELL_SPECIFIC_INCLUDE = {
    "JY": {"ENST00000383126.8": "HLA-DRB3", "ENST00000457451.7": "HLA-DRB4"},
}

CELL_CONFIG = {
    "HCC1937": {
        "factors_table": "uniqueness_correction_L75bp.tsv",
        "conditions": {
            "74_75": ["SRR31477874", "SRR31477875"],   # SRR925709: older batch, not mixed
        },
    },
    "HCT116": {
        "factors_table": "uniqueness_correction_L100bp.tsv",
        "conditions": {
            "41_42": ["SRR3644141", "SRR3644142"],
        },
    },
    "JY": {
        "factors_table": "uniqueness_correction_L100bp.tsv",
        "conditions": {
            "364065": ["SRR364065"],
        },
    },
}

# PIPELINE HOOK: cells/conditions/factor tables from samplesheet
CELL_CONFIG = _CFG.get("CELL_CONFIG", CELL_CONFIG)
CELL_SPECIFIC_INCLUDE = _CFG.get("CELL_SPECIFIC_INCLUDE", CELL_SPECIFIC_INCLUDE)   # PIPELINE HOOK: DR52/DR53 per group
HLA_TOTAL_SOURCE = _CFG.get("hla_total_source", "basic_mane_fc")   # original = basic_mane_fc

OUT_COLS = ['gene_id', 'length', 'counts', 'counts_corrected', 'TPM', 'TPM_corrected',
            'f_unique', 'exclude', 'rule', 'shared']


def load_excluded_set():
    df = pd.read_csv(EXCLUDED_FILE, sep='\t')
    flag = df['excluded'].astype(str).str.upper()
    if not flag.isin(['TRUE', 'FALSE']).all():
        raise SystemExit(f"{EXCLUDED_FILE.name}: `excluded` must be TRUE/FALSE")
    excluded = set(df.loc[flag == 'TRUE', 'transcript_id'])
    included = set(df.loc[flag == 'FALSE', 'transcript_id'])
    if included != set(INCLUDED_ALL):
        raise SystemExit(f"{EXCLUDED_FILE.name} excluded=FALSE rows differ from INCLUDED_ALL: "
                         f"only in table {sorted(included - set(INCLUDED_ALL))}, "
                         f"only in script {sorted(set(INCLUDED_ALL) - included)}")
    dr = set(df.loc[df['DR52_DR53_haplotype'].astype(str).str.upper() == 'TRUE', 'transcript_id'])
    for cell, ids in CELL_SPECIFIC_INCLUDE.items():
        if not set(ids) <= excluded:
            raise SystemExit(f"CELL_SPECIFIC_INCLUDE[{cell}] must list excluded=TRUE transcripts only")
    if not all(set(ids) <= dr for ids in CELL_SPECIFIC_INCLUDE.values()):   # PIPELINE HOOK (was: == JY list)
        raise SystemExit(f"DR52_DR53_haplotype=TRUE rows {sorted(dr)} differ from the JY inclusion list")
    if not set(FOLD_INTO_PRIMARY) <= excluded:
        raise SystemExit("every FOLD_INTO_PRIMARY alt copy must be excluded=TRUE")
    return excluded


def load_pseudogene_set():
    """Transcripts flagged pseudogenes=TRUE in the low-uniqueness annotation table."""
    t = pd.read_csv(LOW_UF_FILE, sep='\t', usecols=['transcript_id', 'pseudogenes'])
    return set(t.loc[t['pseudogenes'].astype(str).str.upper() == 'TRUE', 'transcript_id'])


def load_factors(cell):
    t = pd.read_csv(FACTORS_DIR / CELL_CONFIG[cell]["factors_table"], sep='\t')
    return t.set_index('transcript_id')['uniqueness_factor']


def read_fc(cell, srr, kind):
    """featureCounts table -> DataFrame(length, n) indexed by Geneid."""
    path = BASE / cell / "raw_counts" / f"{srr}_{cell}_{kind}.counts.csv"
    if not path.exists():
        raise SystemExit(f"missing {path}")
    df = pd.read_csv(path, sep='\t', comment='#')
    df = df.rename(columns={'Geneid': 'gene_id', 'Length': 'length', df.columns[-1]: 'n'})
    return df.set_index('gene_id')[['length', 'n']]


def canonical(index):
    """HLA allele IDs (ENST..._A0201) -> canonical MANE ID; others unchanged."""
    return index.str.split('_').str[0]


def correct_sample(cell, srr, factors, excluded_ids, pseudogenes):
    pfc = read_fc(cell, srr, "MANEpersUTR_fc")
    pfcm = read_fc(cell, srr, "MANEpersUTR_fcM")
    if HLA_TOTAL_SOURCE == "basic_mane_fc":                   # original behaviour
        basic = read_fc(cell, srr, "mane_fc")
    else:                                                     # PIPELINE HOOK: one alignment, HLA total = sum of allele fcM
        basic = pd.DataFrame({"length": pfcm["length"].groupby(canonical(pfcm.index)).first(),
                              "n": pfcm["n"].groupby(canonical(pfcm.index)).sum()})

    fc = pfc['n'].groupby(canonical(pfc.index)).sum()
    fcm = pfcm['n'].groupby(canonical(pfcm.index)).sum()
    if set(fcm.index) != set(basic.index) or set(fc.index) != set(basic.index):
        raise SystemExit(f"{srr}: MANEpersUTR and basic MANE transcript sets differ")

    df = basic[['length']].copy()
    df['counts'] = fcm.reindex(df.index)                       # rule 1
    for tid in HLA:                                              # rule 5: HLA total from basic MANE
        df.loc[tid, 'counts'] = basic.loc[tid, 'n']
    fcm_a, fc_a = fcm.reindex(df.index), fc.reindex(df.index)   # align on transcript ID
    df['shared'] = (1 - fc_a / fcm_a).where(fcm_a >= MIN_FCM_FOR_SHARED)
    # HLA: in MANEpersUTR the sample's two alleles compete with each other, so fc/fcM
    # does not measure competition with other genes -> treated as "not measurable"
    df.loc[list(HLA), 'shared'] = np.nan

    # raw TPM from the uncorrected counts over all transcripts
    rpk = df['counts'] / (df['length'] / 1000.0)
    df['TPM'] = rpk / rpk.sum() * 1_000_000

    # rule 2: fold alt duplicates into their primary copy
    work = df['counts'].copy()
    for alt, prim in FOLD_INTO_PRIMARY.items():
        work[prim] += work[alt]

    # rule 6: exclusions
    df['exclude'] = df.index.isin(excluded_ids)
    df['rule'] = ''
    df.loc[df.index.isin(excluded_ids), 'rule'] = 'excluded: non-primary'
    df.loc[df.index.isin(list(FOLD_INTO_PRIMARY)), 'rule'] = 'excluded: alt copy folded into primary'
    df.loc[df.index.isin(list(PAR_GENES)), ['exclude', 'rule']] = [True, 'excluded: PAR']
    nc = set(pd.read_csv(NONCODING_FILE, sep='\t')['transcript_id'])
    df.loc[df.index.isin(nc), ['exclude', 'rule']] = [True, 'excluded: non-coding']

    # rules 3-5: f_unique
    uf = factors.reindex(df.index)
    df['f_unique'] = 1.0
    keep = ~df['exclude']
    df.loc[keep, 'rule'] = 'f=1'
    rule3 = (keep & (uf < UF_MAX) & df.index.isin(pseudogenes)
             & (df['shared'].isna() | (df['shared'] < SHARED_MAX)))
    df.loc[rule3, 'f_unique'] = uf[rule3].where(uf[rule3] > 0, ZERO_FLOOR)
    df.loc[rule3, 'rule'] = 'pseudogene inflation: f=UF'
    mane_comp = keep & (uf < UF_MAX) & df.index.isin(pseudogenes) & (df['shared'] >= SHARED_MAX)
    df.loc[mane_comp, 'rule'] = 'f=1: MANE competition (fcM>fc)'
    df.loc[df.index.isin(list(FOLD_INTO_PRIMARY.values())), 'rule'] += ' (+alt fcM)'
    df.loc[df['exclude'], 'f_unique'] = np.nan

    df['counts_corrected'] = np.where(df['exclude'], np.nan, work * df['f_unique'])
    rpk_c = df['counts_corrected'] / (df['length'] / 1000.0)
    df['TPM_corrected'] = rpk_c / rpk_c.sum(skipna=True) * 1_000_000

    # rule 5: split canonical HLA rows into the sample's alleles (heterozygous only)
    df = df.reset_index()
    rows = []
    for tid, gene in HLA.items():
        alleles = pfc.loc[pfc.index.str.startswith(tid + '_'), 'n']
        if alleles.empty:
            continue                                             # homozygous: canonical row stays
        if alleles.sum() == 0:
            raise SystemExit(f"{srr}: no unique {gene} allele reads to derive a ratio")
        ratio = alleles / alleles.sum()
        base = df[df['gene_id'] == tid].iloc[0]
        for aid, r in ratio.items():
            row = base.copy()
            row['gene_id'] = aid
            for c in ['counts', 'counts_corrected', 'TPM', 'TPM_corrected']:
                row[c] = base[c] * r
            row['rule'] = f"{base['rule']}; allele share {r:.3f} (unique MANEpersUTR reads)"
            rows.append(row)
        df = df[df['gene_id'] != tid]
    if rows:
        df = pd.concat([df, pd.DataFrame(rows)], ignore_index=True)

    print(f"    [{srr}] excluded={int(df['exclude'].sum())}  "
          f"rule3 (f=UF)={int(rule3.sum())}  MANE-competition f=1={int(mane_comp.sum())}  "
          f"TPM_corrected_sum={df['TPM_corrected'].sum(skipna=True):,.0f}")
    return df[OUT_COLS]


def average_condition(sample_dfs):
    """Mean counts/counts_corrected/TPM/TPM_corrected across replicate samples,
    keyed by gene_id. Excluded transcripts stay excluded (NaN) in the mean."""
    merged = None
    for i, df in enumerate(sample_dfs):
        d = df[['gene_id', 'length', 'counts', 'counts_corrected', 'TPM', 'TPM_corrected']].copy()
        d = d.rename(columns={c: f"{c}__{i}" for c in
                               ['counts', 'counts_corrected', 'TPM', 'TPM_corrected']})
        merged = d if merged is None else merged.merge(d, on=['gene_id', 'length'], how='outer')

    n = len(sample_dfs)
    out = merged[['gene_id', 'length']].copy()
    for col in ['counts', 'counts_corrected', 'TPM', 'TPM_corrected']:
        cols = [f"{col}__{i}" for i in range(n)]
        out[f"mean_{col}"] = merged[cols].mean(axis=1, skipna=False)
    return out.sort_values('gene_id')


def main():
    excluded_ids = load_excluded_set()
    pseudogenes = load_pseudogene_set()
    print(f"Excluded non-primary transcripts: {len(excluded_ids)}  |  PAR: {len(PAR_GENES)}  |  "
          f"pseudogene-flagged (low-UF table): {len(pseudogenes)}")

    for cell, cfg in CELL_CONFIG.items():
        print(f"\n{'=' * 60}\n{cell}\n{'=' * 60}")
        cell_excluded = excluded_ids - set(CELL_SPECIFIC_INCLUDE.get(cell, {}))
        factors = load_factors(cell)

        outdir = BASE / cell / "corrected_counts"
        outdir.mkdir(exist_ok=True)

        for condition, srrs in cfg["conditions"].items():
            print(f"\n  Condition: {condition}  (samples: {', '.join(srrs)})")
            sample_dfs = []
            for srr in srrs:
                df = correct_sample(cell, srr, factors, cell_excluded, pseudogenes)
                df.to_csv(outdir / f"{srr}_mane_fc_corrected.tsv",
                          sep='\t', index=False, float_format='%.6f')
                sample_dfs.append(df)

            mean_df = average_condition(sample_dfs)
            mean_path = outdir / f"{cell}_corrected_mean_{condition}_mane_fc_corrected.tsv"
            mean_df.to_csv(mean_path, sep='\t', index=False, float_format='%.6f')
            print(f"    -> {mean_path}  ({len(mean_df):,} transcripts, "
                  f"mean_TPM_corrected_sum={mean_df['mean_TPM_corrected'].sum(skipna=True):,.0f})")

    print("\nDone.")


if __name__ == '__main__':
    main()
