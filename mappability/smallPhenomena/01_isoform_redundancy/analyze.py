#!/usr/bin/env python3
"""
01_isoform_redundancy/analyze.py  [v2]

Central question: when a gene has multiple annotated CDS isoforms or belongs to a
paralog family, how much signal can short-read tools actually resolve — and does the
UTR-CDS-UTR transcript structure help?

Analyses
--------
A  isoform_table.tsv          — n_isoforms, lengths per gene (run once)
B  isoform_groups.tsv         — + pseudogene count + functional groups (run once)
C  isoform_kmers_cds.tsv      — unique 31-mer fraction + max unique stretch (bp), CDS (run once)
   isoform_kmers_mane.tsv     — unique 31-mer fraction + max unique stretch (bp), MANE (run once)
D  runs/{label}/vs_{index}/   — all 22 kallisto conditions × 2 indices + v8 proper
E  isoform_homology.tsv       — MAFFT pairwise CDS isoform identity (run once)
F  utr_fail_stats.tsv         — χ², entropy, CV for UTR-unrescued genes
G  paralog_homology_mane.tsv  — MAFFT on paralog families (MANE canonical seqs)
   paralog_homology_isoforms.tsv  — same but all GENCODE isoforms (≤20/member)
   id_verification.tsv        — cross-source ID consistency check
H  gene_level_recommendation.tsv   — isoform-redundant, collapse to gene level
   cluster_level_recommendation.tsv — paralog-ambiguous, work at cluster level

Run registry (Section D)
------------------------
22 read conditions × 2 kallisto indices (vs_dedup_cds, vs_MANE):
  dedup_cds_reads  {PE,SE} × {75,100,150,200} bp  — 8 conditions
  MANE_RNA_reads   {PE,SE} × {75,100,150,200} bp  — 8 conditions
  MANE_exonflank   {PE,SE} × {75,100,150}    bp   — 6 conditions

v8 proper-pair (3 STAR-UF conditions, ENST-keyed, MANE RNA reads 150bp PE):
  v8_proper_genomic, v8_proper_dedup_cds, v8_proper_MANE
"""

import csv
import gzip
import re
import subprocess
import sys
from collections import Counter, defaultdict
from math import log2
from pathlib import Path
from statistics import mean, median, stdev

# scipy is only required for Section F (chi-squared test)
try:
    from scipy.stats import chisquare as scipy_chisquare
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print('[warn] scipy not found — chi-squared test will be skipped (install scipy)')

HERE   = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from lib_check import check_library, warn_if_mixing_runs

ROOT    = HERE.parent.parent.parent          # mappability-correction/
MAPDIR  = ROOT / 'mappability'
REFDIR  = ROOT / 'Ref'
MANEDIR = Path(__file__).resolve().parents[4] / 'MANE'
GENDIR  = Path(__file__).resolve().parents[4] / 'GENCODE'

CDS_FA      = REFDIR / 'GCF_000001405.40_GRCh38.p14_cds_from_genomic.dedup.fna'
MANE_FA     = REFDIR / 'MANE.GRCh38.v1.5.ensembl_rna.fna'
GENCODE_FA  = GENDIR / 'gencode.v49.pc_transcripts.fa'   # full UTR-CDS-UTR all isoforms
CONV_TSV    = MAPDIR / 'results' / 'id_conversion_table.tsv'
GENCODE_GTF = GENDIR / 'gencode.v49.annotation.gtf.gz'

GROUP_FILES = {
    'kinase':        MANEDIR / 'kinase_enst.txt',
    'tf':            MANEDIR / 'tf_enst.txt',
    'cell_cycle':    MANEDIR / 'cell_cycle_enst.txt',
    'extracellular': MANEDIR / 'combined_extracellular_enst.txt',
    'blacklist':     MANEDIR / 'functional_blacklist_enst.txt',
}

KMER_K          = 31
MAFFT_MAX_LEN   = 4_000   # nt — CDS isoform homology cap (Section E)
MAFFT_MAX_GENES = 500     # gene cap for Section E
UTR_FAIL_THRESH = 0.90    # frac_unique threshold below which UTR rescue is insufficient
PARALOG_MAX     = 99      # max family members for MAFFT (Section G)
PARALOG_ISOFORM_MAX = 20  # max isoforms per family member for option B

# ── Run registry ──────────────────────────────────────────────────────────────
# Columns: (label, mode, read_len_bp, read_source, star_uf_fname_or_None, kallisto_dir)
#
# read_source values:
#   'dedup_cds'    — reads simulated from deduplicated CDS sequences (no UTRs)
#   'MANE_RNA'     — reads from full spliced MANE transcripts (5'UTR+CDS+3'UTR)
#   'MANE_exonflank' — reads spanning exon–intron boundaries (±50 bp flank)
#
# star_uf_fname: CDS-level selfmap file in results/ — only available for
#   dedup_cds PE/SE 100 bp (the only CDS selfmap STAR runs executed so far).
#   When None, STAR columns are left empty in isoform_mappability.tsv.
#
# kallisto_dir: basename under kallisto_results/ — each dir contains
#   kallisto_vs_dedup_cds/abundance.tsv  (dedup CDS index, NP_ target IDs)
#   kallisto_vs_MANE/abundance.tsv       (MANE index, ENST target IDs)

_RUNS = [
    # ── dedup CDS reads ───────────────────────────────────────────────────────
    # Reads come from CDS sequences only (no UTRs). The dedup-CDS selfmap is the
    # floor: how much of each isoform is distinguishable by CDS k-mers alone.
    # vs_MANE shows how much UTR k-mers in the MANE index recover vs the CDS floor.
    ('dedup_cds_reads_PE_L75bp',  'PE',  75,  'dedup_cds', None,
     'mappability_dedup_cds_reads_PE_L75bp'),
    ('dedup_cds_reads_PE_L100bp', 'PE', 100,  'dedup_cds',
     'cds_uniqueness_factors_selfmap_PE.tsv',
     'mappability_dedup_cds_reads_PE'),
    ('dedup_cds_reads_PE_L150bp', 'PE', 150,  'dedup_cds', None,
     'mappability_dedup_cds_reads_PE_L150bp'),
    ('dedup_cds_reads_PE_L200bp', 'PE', 200,  'dedup_cds', None,
     'mappability_dedup_cds_reads_PE_L200bp'),
    ('dedup_cds_reads_SE_L75bp',  'SE',  75,  'dedup_cds', None,
     'mappability_dedup_cds_reads_SE_L75bp'),
    ('dedup_cds_reads_SE_L100bp', 'SE', 100,  'dedup_cds',
     'cds_uniqueness_factors_selfmap_SE.tsv',
     'mappability_dedup_cds_reads_SE'),
    ('dedup_cds_reads_SE_L150bp', 'SE', 150,  'dedup_cds', None,
     'mappability_dedup_cds_reads_SE_L150bp'),
    ('dedup_cds_reads_SE_L200bp', 'SE', 200,  'dedup_cds', None,
     'mappability_dedup_cds_reads_SE_L200bp'),

    # ── MANE RNA reads ────────────────────────────────────────────────────────
    # Reads from full spliced MANE transcripts (5'UTR + CDS + 3'UTR). This is
    # the closest proxy to real poly-A RNA-seq.
    # vs_dedup_cds shows how many full-transcript reads survive CDS-only quantification.
    # vs_MANE shows gene-level mappability with UTR disambiguation.
    ('MANE_RNA_reads_PE_L75bp',  'PE',  75,  'MANE_RNA', None,
     'mappability_MANE_RNA_PE_L75bp'),
    ('MANE_RNA_reads_PE_L100bp', 'PE', 100,  'MANE_RNA', None,
     'mappability_MANE_RNA_PE'),
    ('MANE_RNA_reads_PE_L150bp', 'PE', 150,  'MANE_RNA', None,
     'mappability_MANE_RNA_PE_L150bp'),
    ('MANE_RNA_reads_PE_L200bp', 'PE', 200,  'MANE_RNA', None,
     'mappability_MANE_RNA_PE_L200bp'),
    ('MANE_RNA_reads_SE_L75bp',  'SE',  75,  'MANE_RNA', None,
     'mappability_MANE_RNA_SE_L75bp'),
    ('MANE_RNA_reads_SE_L100bp', 'SE', 100,  'MANE_RNA', None,
     'mappability_MANE_RNA_SE'),
    ('MANE_RNA_reads_SE_L150bp', 'SE', 150,  'MANE_RNA', None,
     'mappability_MANE_RNA_SE_L150bp'),
    ('MANE_RNA_reads_SE_L200bp', 'SE', 200,  'MANE_RNA', None,
     'mappability_MANE_RNA_SE_L200bp'),

    # ── MANE exon-flank reads ─────────────────────────────────────────────────
    # Reads spanning exon–intron boundaries (±50 bp). Approximates total RNA-seq:
    # exon-only reads pseudoalign normally; junction-spanning reads fail the spliced
    # index. The failure rate is a lower-bound on total-RNA batch effects.
    ('MANE_exonflank_reads_PE_L75bp',  'PE',  75,  'MANE_exonflank', None,
     'mappability_MANE_exonflank_PE_L75bp'),
    ('MANE_exonflank_reads_PE_L100bp', 'PE', 100,  'MANE_exonflank', None,
     'mappability_MANE_exonflank_PE'),
    ('MANE_exonflank_reads_PE_L150bp', 'PE', 150,  'MANE_exonflank', None,
     'mappability_MANE_exonflank_PE_L150bp'),
    ('MANE_exonflank_reads_SE_L75bp',  'SE',  75,  'MANE_exonflank', None,
     'mappability_MANE_exonflank_SE_L75bp'),
    ('MANE_exonflank_reads_SE_L100bp', 'SE', 100,  'MANE_exonflank', None,
     'mappability_MANE_exonflank_SE'),
    ('MANE_exonflank_reads_SE_L150bp', 'SE', 150,  'MANE_exonflank', None,
     'mappability_MANE_exonflank_SE_L150bp'),
]

# v8 proper-pair: MANE RNA reads (150 bp PE) with correct 300 bp fragment geometry.
# R1 and R2 come from opposite ends of the fragment (unlike v1–v7 where both reads
# covered the same window). These use ENST IDs (not NP_ CDS IDs) and are handled
# by analyse_v8_proper() rather than analyse_run().
_V8_PROPER = [
    ('v8_proper_genomic',
     'transcript_uniqueness_factors_genomic_RNA_PE_L150bp_proper.tsv',
     'Full genome GRCh38.p14'),
    ('v8_proper_dedup_cds',
     'transcript_uniqueness_factors_dedup_cds_RNA_PE_L150bp_proper.tsv',
     'Deduplicated CDS reference'),
    ('v8_proper_MANE',
     'transcript_uniqueness_factors_MANE_RNA_PE_L150bp_proper.tsv',
     'MANE Select transcriptome'),
]


# ══════════════════════════════════════════════════════════════════════════════
# FASTA parsers
# ══════════════════════════════════════════════════════════════════════════════

def parse_cds_fasta(path):
    """
    Parse deduplicated CDS FASTA.
    Header: >lcl|... [gene=SYMBOL] [protein_id=NP_...]
    Returns: gene2ids {symbol: [cds_id]}, meta {cds_id: (gene, prot, loc, len)},
             seqs {cds_id: sequence}
    """
    gene2ids = defaultdict(list)
    meta, seqs = {}, {}
    cur_id, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if cur_id:
                    seqs[cur_id] = ''.join(buf)
                    meta[cur_id] = meta[cur_id][:3] + (len(seqs[cur_id]),)
                cur_id = line[1:].split()[0]
                gm = re.search(r'\[gene=([^\]]+)\]', line)
                pm = re.search(r'\[protein_id=([^\]]+)\]', line)
                lm = re.search(r'\[location=([^\]]+)\]', line)
                gene = gm.group(1) if gm else ''
                meta[cur_id] = (gene, pm.group(1) if pm else '',
                                lm.group(1) if lm else '', 0)
                if gene:
                    gene2ids[gene].append(cur_id)
                buf = []
            else:
                buf.append(line.upper())
    if cur_id:
        seqs[cur_id] = ''.join(buf)
        meta[cur_id] = meta[cur_id][:3] + (len(seqs[cur_id]),)
    return gene2ids, meta, seqs


def parse_mane_fasta(path):
    """
    Parse MANE Select RNA FASTA (full UTR-CDS-UTR sequences).
    Header: >ENST... cdna ... gene_symbol:SYMBOL
    Returns: symbol2enst {symbol: enst_id}, enst2seq {enst_id: seq}
    """
    symbol2enst, enst2seq = {}, {}
    cur_id, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if cur_id:
                    enst2seq[cur_id] = ''.join(buf)
                parts = line[1:].split()
                cur_id = parts[0]
                sym_m = re.search(r'gene_symbol:(\S+)', line)
                if sym_m:
                    symbol2enst[sym_m.group(1)] = cur_id
                buf = []
            else:
                buf.append(line.upper())
    if cur_id:
        enst2seq[cur_id] = ''.join(buf)
    return symbol2enst, enst2seq


def parse_gencode_pc_fasta(path):
    """
    Parse GENCODE v49 protein-coding transcripts FASTA (all isoforms, UTR-CDS-UTR).
    Header: >ENST|ENSG|...|transcript_name|SYMBOL|length|UTR5:...|CDS:...|UTR3:...|
    Gene symbol is pipe-delimited field index 5 (0-based).
    Returns: symbol2ensts {symbol: [enst_id]}, enst2seq {enst_id: seq},
             enst2symbol {enst_id: symbol}
    """
    symbol2ensts = defaultdict(list)
    enst2seq, enst2symbol = {}, {}
    cur_id, cur_sym, buf = None, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if cur_id:
                    enst2seq[cur_id] = ''.join(buf)
                fields = line[1:].split('|')
                cur_id = fields[0] if fields else None
                # Field 5 is the gene symbol (e.g. OR4F5)
                cur_sym = fields[5] if len(fields) > 5 else ''
                if cur_id and cur_sym:
                    symbol2ensts[cur_sym].append(cur_id)
                    enst2symbol[cur_id] = cur_sym
                buf = []
            else:
                buf.append(line.upper())
    if cur_id:
        enst2seq[cur_id] = ''.join(buf)
    return symbol2ensts, enst2seq, enst2symbol


# ══════════════════════════════════════════════════════════════════════════════
# A. Isoform table
# ══════════════════════════════════════════════════════════════════════════════

def make_isoform_table(gene2ids, meta):
    rows = []
    for gene in sorted(gene2ids):
        ids = gene2ids[gene]
        lengths = [meta[i][3] for i in ids]
        rows.append({
            'gene_name':          gene,
            'n_isoforms':         len(ids),
            'total_length_nt':    sum(lengths),
            'min_length_nt':      min(lengths),
            'max_length_nt':      max(lengths),
            'protein_ids':        ';'.join(meta[i][1] for i in ids),
            'isoform_ids':        ';'.join(ids),
            'isoform_lengths_nt': ';'.join(str(l) for l in lengths),
        })
    return rows


def write_isoform_table(rows):
    fields = ['gene_name', 'n_isoforms', 'total_length_nt', 'min_length_nt',
              'max_length_nt', 'protein_ids', 'isoform_ids', 'isoform_lengths_nt']
    out = HERE / 'isoform_table.tsv'
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    print(f'  Saved {out.name}  ({len(rows):,} genes)')
    return rows


# ══════════════════════════════════════════════════════════════════════════════
# B. Functional groups
# ══════════════════════════════════════════════════════════════════════════════

def load_id_conversion(path):
    """Returns np2sym {NP_: symbol} and enst2sym {ENST: symbol}."""
    np2sym, enst2sym = {}, {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            sym  = row['gene_symbol']
            enst = row['ensembl_transcript_id']
            np_  = row['refseq_protein_id']
            if enst: enst2sym[enst] = sym
            if np_:  np2sym[np_]   = sym
    return np2sym, enst2sym


def load_groups(enst2sym):
    """Returns {symbol: set of group labels}."""
    sym2groups = defaultdict(set)
    for group, path in GROUP_FILES.items():
        if not path.exists():
            print(f'  [warn] group file missing: {path.name}')
            continue
        with open(path) as fh:
            for line in fh:
                enst = line.strip().split('\t')[0]
                if enst in enst2sym:
                    sym2groups[enst2sym[enst]].add(group)
    return sym2groups


def count_pseudogenes_gencode(coding_genes):
    """Assign pseudogenes to parent coding genes by gene_name prefix matching."""
    print('  Parsing GENCODE v49 for pseudogene counts ...')
    pseudo_biotypes = {
        'processed_pseudogene', 'unprocessed_pseudogene',
        'transcribed_processed_pseudogene', 'transcribed_unprocessed_pseudogene',
        'transcribed_unitary_pseudogene', 'unitary_pseudogene',
        'rRNA_pseudogene', 'IG_V_pseudogene', 'TR_V_pseudogene', 'TR_J_pseudogene',
    }
    pseudo_names = []
    opener = gzip.open if str(GENCODE_GTF).endswith('.gz') else open
    with opener(GENCODE_GTF, 'rt') as fh:
        for line in fh:
            if line.startswith('#'): continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'gene': continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', f[8]))
            if attrs.get('gene_type') in pseudo_biotypes:
                pseudo_names.append(attrs.get('gene_name', ''))
    print(f'  {len(pseudo_names):,} pseudogene gene entries found')

    coding_sorted = sorted(coding_genes, key=len, reverse=True)
    pseudo_counts = Counter()
    for pname in pseudo_names:
        if not pname: continue
        pu = pname.upper()
        for sym in coding_sorted:
            su = sym.upper()
            if pu.startswith(su) and len(pu) > len(su):
                nxt = pu[len(su)]
                if nxt in ('P', '_', '-'):
                    pseudo_counts[sym] += 1
                    break
    matched = sum(pseudo_counts.values())
    print(f'  {matched:,} pseudogenes matched to {len(pseudo_counts):,} coding genes '
          f'({100*matched/max(len(pseudo_names),1):.1f}%)')
    return pseudo_counts


def write_isoform_groups(iso_rows, pseudo_counts, sym2groups, np2sym, meta, gene2ids):
    def isoform_bin(n):
        if n == 1:  return '1'
        if n <= 3:  return '2-3'
        if n <= 5:  return '4-5'
        if n <= 10: return '6-10'
        return '10+'

    rows = []
    for r in iso_rows:
        sym = r['gene_name']
        grp = sym2groups.get(sym, set())
        rows.append({
            'gene_name':        sym,
            'n_isoforms':       r['n_isoforms'],
            'pseudogene_count': pseudo_counts.get(sym, 0),
            'is_kinase':        int('kinase' in grp),
            'is_tf':            int('tf' in grp),
            'is_cell_cycle':    int('cell_cycle' in grp),
            'is_extracellular': int('extracellular' in grp),
            'in_blacklist':     int('blacklist' in grp),
            'min_isoform_len_nt': r['min_length_nt'],
            'max_isoform_len_nt': r['max_length_nt'],
            'n_isoforms_bin':   isoform_bin(r['n_isoforms']),
        })

    fields = ['gene_name', 'n_isoforms', 'pseudogene_count',
              'is_kinase', 'is_tf', 'is_cell_cycle', 'is_extracellular', 'in_blacklist',
              'min_isoform_len_nt', 'max_isoform_len_nt', 'n_isoforms_bin']
    out = HERE / 'isoform_groups.tsv'
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    print(f'  Saved {out.name}  ({len(rows):,} genes)')

    print('\n  ── n_isoforms vs pseudogene_count (by isoform bin) ──')
    print(f'  {"bin":<8} {"n_genes":>8} {"mean_pseudo":>12} {"mean_n_iso":>11}')
    bins = defaultdict(lambda: {'iso': [], 'ps': []})
    for r in rows:
        b = r['n_isoforms_bin']
        bins[b]['iso'].append(r['n_isoforms'])
        bins[b]['ps'].append(r['pseudogene_count'])
    for b in ['1', '2-3', '4-5', '6-10', '10+']:
        if b not in bins: continue
        d = bins[b]
        print(f'  {b:<8} {len(d["iso"]):>8,} '
              f'{mean(d["ps"]):>12.2f} {mean(d["iso"]):>11.2f}')

    single = [r for r in rows if r['n_isoforms'] == 1]
    n_tot  = len(rows)
    n_s    = len(single)
    print(f'\n  ── Functional group enrichment in single-isoform genes ──')
    print(f'  Single-isoform: {n_s:,} / {n_tot:,} ({100*n_s/n_tot:.1f}%)')
    for g in ('kinase', 'tf', 'cell_cycle', 'extracellular', 'blacklist'):
        col = f'is_{g}' if g != 'blacklist' else 'in_blacklist'
        n_g  = sum(1 for r in rows if r[col])
        n_gs = sum(1 for r in single if r[col])
        if n_g == 0: continue
        print(f'  {g:<15} {n_gs:>5}/{n_g:<5} '
              f'({100*n_gs/n_g:.1f}% of group are single-isoform; '
              f'overall rate {100*n_s/n_tot:.1f}%)')
    return rows


# ══════════════════════════════════════════════════════════════════════════════
# C. 31-mer unique fraction + unique stretch analysis
# ══════════════════════════════════════════════════════════════════════════════

MIN_UNIQUE_STRETCH_BP = 50   # target: unique stretch spans at least this many bp
# A stretch of N consecutive unique 31-mer positions spans N + k - 1 bp.
# To reach MIN_UNIQUE_STRETCH_BP = 50: need N ≥ 50 - 31 + 1 = 20 consecutive positions.
# One unique 31-mer (N=1) spans 31 bp — the minimum reportable case.
# Rationale: a 50 bp unique region can accommodate short indels or mismatches
# that would shift a read's k-mer positions by a few bp and still land on
# at least one unique 31-mer.


def kmers(seq, k=KMER_K):
    return {seq[i:i+k] for i in range(len(seq) - k + 1) if 'N' not in seq[i:i+k]}


def kmers_positional(seq, k=KMER_K):
    """Return list of (position, kmer) for all valid (N-free) k-mers."""
    return [(i, seq[i:i+k]) for i in range(len(seq) - k + 1)
            if 'N' not in seq[i:i+k]]


def max_unique_stretch(seq, unique_set, k=KMER_K):
    """
    Find the longest contiguous run of positions i where seq[i:i+k] is in unique_set.
    Positions are consecutive integers (i, i+1, i+2, ...) — each step adds 1 bp.
    The span in bp of a run of N positions = N + k - 1.

    Returns (max_stretch_bp, max_stretch_positions, best_start) where:
      max_stretch_bp        : span in bp of the longest unique run (31 if N=1)
      max_stretch_positions : number of consecutive positions in the best run
      best_start            : 0-based start position of that run in the sequence
    """
    best_len, best_start = 0, 0
    cur_len,  cur_start  = 0, 0
    for i, km in kmers_positional(seq, k):
        if km in unique_set:
            if cur_len == 0:
                cur_start = i
            cur_len += 1
        else:
            if cur_len > best_len:
                best_len, best_start = cur_len, cur_start
            cur_len = 0
    if cur_len > best_len:
        best_len, best_start = cur_len, cur_start
    stretch_bp = best_len + k - 1 if best_len > 0 else 0
    return stretch_bp, best_len, best_start


def kmer_unique_fraction_within_gene(gene2ids, seqs, k=KMER_K):
    """
    Within-gene k-mer uniqueness for multi-isoform CDS genes.
    'unique' here means the k-mer is not shared with ANY other isoform of the same gene.
    Also reports the longest contiguous run of unique k-mer positions per isoform.
    """
    rows = []
    for gene in sorted(gene2ids):
        ids = gene2ids[gene]
        if len(ids) < 2: continue
        iso_seqs = {i: seqs[i] for i in ids if i in seqs and len(seqs[i]) >= k}
        if len(iso_seqs) < 2: continue
        iso_kmers = {i: kmers(s, k) for i, s in iso_seqs.items()}
        freq = Counter(km for km_set in iso_kmers.values() for km in km_set)
        for iso_id, km_set in iso_kmers.items():
            total      = len(km_set)
            unique_set = {km for km in km_set if freq[km] == 1}
            unique     = len(unique_set)
            stretch_bp, stretch_pos, stretch_start = max_unique_stretch(
                iso_seqs[iso_id], unique_set, k)
            rows.append({
                'gene_name':             gene,
                'n_isoforms':            len(ids),
                'isoform_id':            iso_id,
                'total_kmers':           total,
                'unique_kmers':          unique,
                'shared_kmers':          total - unique,
                'frac_unique':           round(unique/total, 4) if total else 0.0,
                'frac_shared':           round((total-unique)/total, 4) if total else 1.0,
                'max_unique_stretch_bp': stretch_bp,
                'max_stretch_positions': stretch_pos,
                'stretch_start_nt':      stretch_start,
                'stretch_reaches_50bp':  int(stretch_bp >= MIN_UNIQUE_STRETCH_BP),
                'reference':             'CDS',
                'k':                     k,
            })
    return rows


def kmer_unique_fraction_mane(symbol2enst, enst2seq, k=KMER_K):
    """
    Cross-gene k-mer uniqueness over all MANE transcripts (UTR + CDS + UTR).
    'unique' means the k-mer does not appear in any OTHER MANE transcript.
    Two-pass: build global frequency, then score each transcript.
    Also reports the longest contiguous run of unique k-mer positions.
    """
    print(f'  Building global {k}-mer frequency across {len(enst2seq):,} MANE transcripts ...')
    freq = Counter()
    for seq in enst2seq.values():
        freq.update(kmers(seq, k))

    enst2sym = {v: kk for kk, v in symbol2enst.items()}
    rows = []
    for enst, seq in sorted(enst2seq.items()):
        if len(seq) < k: continue
        all_km     = kmers(seq, k)
        total      = len(all_km)
        unique_set = {m for m in all_km if freq[m] == 1}
        unique     = len(unique_set)
        stretch_bp, stretch_pos, stretch_start = max_unique_stretch(seq, unique_set, k)
        rows.append({
            'gene_name':             enst2sym.get(enst, ''),
            'enst_id':               enst,
            'seq_len_nt':            len(seq),
            'total_kmers':           total,
            'unique_kmers':          unique,
            'shared_kmers':          total - unique,
            'frac_unique':           round(unique/total, 4) if total else 0.0,
            'frac_shared':           round((total-unique)/total, 4) if total else 1.0,
            'max_unique_stretch_bp': stretch_bp,
            'max_stretch_positions': stretch_pos,
            'stretch_start_nt':      stretch_start,
            'stretch_reaches_50bp':  int(stretch_bp >= MIN_UNIQUE_STRETCH_BP),
            'reference':             'MANE_RNA',
            'k':                     k,
        })
    return rows


def write_kmer_table(rows, out_path):
    if not rows: return
    fields = list(rows[0].keys())
    with open(out_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    frac    = [r['frac_unique'] for r in rows]
    n_reach = sum(r['stretch_reaches_50bp'] for r in rows)
    print(f'  Saved {out_path.name}  ({len(rows):,} entries)  '
          f'mean frac_unique={mean(frac):.3f}  median={median(frac):.3f}  '
          f'stretch≥50bp: {n_reach:,}/{len(rows):,} ({100*n_reach/len(rows):.1f}%)')


# ══════════════════════════════════════════════════════════════════════════════
# D. Per-run STAR UF + kallisto TPM (all 22 conditions × 2 indices + v8 proper)
# ══════════════════════════════════════════════════════════════════════════════

def load_star_cds_uf(path):
    """Load CDS-level selfmap UF. Returns {cds_id: uf}."""
    d = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            d[row['cds_id']] = float(row['uniqueness_factor'])
    return d


def load_kallisto(path):
    """
    Load kallisto abundance.tsv.
    Returns {target_id: {'tpm': float, 'est_counts': float}}.
    Both TPM and est_counts are kept: TPM for cross-sample comparison,
    est_counts for the chi-squared test in Section F.
    Note: est_counts are EM-estimated fractions, NOT raw integer counts.
    Do not compare est_counts across runs without TPM-normalisation.
    """
    d = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            d[row['target_id']] = {
                'tpm':        float(row['tpm']),
                'est_counts': float(row['est_counts']),
            }
    return d


def entropy(values):
    """Normalised Shannon entropy [0, 1]. Returns 0 if len(values) < 2."""
    total = sum(values)
    if total == 0 or len(values) < 2: return 0.0
    probs = [v/total for v in values]
    raw   = -sum(p * log2(p) for p in probs if p > 0)
    return raw / log2(len(values))


def _write_run_outputs(outdir, gene2ids, meta, star_uf, kal_data,
                       rlen, mode, index_label):
    """
    Write isoform_mappability.tsv and reads_by_nisoform.tsv to outdir.

    star_uf   : {cds_id: float} or {} if no STAR selfmap for this condition
    kal_data  : {target_id: {'tpm', 'est_counts'}}
    index_label : 'dedup_cds' or 'MANE'

    When index_label == 'MANE', target IDs are ENST — one per gene.
    The n_with_kallisto column will be 1 for every gene; entropy is always 0
    (nothing to distribute), but kallisto_tpm_sum captures gene-level TPM
    with UTR-based disambiguation.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    map_rows = []
    for gene in sorted(gene2ids):
        ids  = gene2ids[gene]
        ufs  = [star_uf[i] for i in ids if i in star_uf]
        # For dedup_cds index: targets are NP_ CDS IDs directly.
        # For MANE index: targets are ENST IDs — match via gene name.
        tpm_vals    = []
        count_vals  = []
        for i in ids:
            if i in kal_data:
                tpm_vals.append(kal_data[i]['tpm'])
                count_vals.append(kal_data[i]['est_counts'])

        map_rows.append({
            'gene_name':            gene,
            'n_isoforms':           len(ids),
            'kallisto_index':       index_label,
            'n_with_star_uf':       len(ufs),
            'star_uf_mean':         round(mean(ufs), 4)   if ufs else '',
            'star_uf_median':       round(median(ufs), 4) if ufs else '',
            'star_uf_min':          round(min(ufs), 4)    if ufs else '',
            'star_uf_max':          round(max(ufs), 4)    if ufs else '',
            'star_n_uf_zero':       sum(1 for u in ufs if u == 0.0) if ufs else '',
            'star_frac_uf_zero':    round(sum(1 for u in ufs if u == 0.0)/len(ufs), 4) if ufs else '',
            'n_with_kallisto':      len(tpm_vals),
            'kallisto_tpm_sum':     round(sum(tpm_vals), 3)        if tpm_vals else '',
            'kallisto_tpm_max':     round(max(tpm_vals), 3)        if tpm_vals else '',
            'kallisto_tpm_entropy': round(entropy(tpm_vals), 4)    if tpm_vals else '',
            'kallisto_tpm_cv':      (round(stdev(tpm_vals)/mean(tpm_vals), 4)
                                     if len(tpm_vals) > 1 and mean(tpm_vals) > 0 else ''),
            'kallisto_est_counts_sum': round(sum(count_vals), 3)   if count_vals else '',
            'read_len_bp':          rlen,
            'mode':                 mode,
        })

    map_fields = [
        'gene_name', 'n_isoforms', 'kallisto_index',
        'n_with_star_uf', 'star_uf_mean', 'star_uf_median', 'star_uf_min', 'star_uf_max',
        'star_n_uf_zero', 'star_frac_uf_zero',
        'n_with_kallisto', 'kallisto_tpm_sum', 'kallisto_tpm_max',
        'kallisto_tpm_entropy', 'kallisto_tpm_cv', 'kallisto_est_counts_sum',
        'read_len_bp', 'mode',
    ]
    with open(outdir / 'isoform_mappability.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=map_fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(map_rows)

    # Bin by isoform count — only meaningful for multi-isoform genes
    bins = defaultdict(lambda: {'ufs': [], 'ents': [], 'fzeros': []})
    for r in map_rows:
        n = min(r['n_isoforms'], 10)
        if r['star_uf_mean']        != '': bins[n]['ufs'].append(float(r['star_uf_mean']))
        if r['kallisto_tpm_entropy'] != '': bins[n]['ents'].append(float(r['kallisto_tpm_entropy']))
        if r['star_frac_uf_zero']   != '': bins[n]['fzeros'].append(float(r['star_frac_uf_zero']))

    bin_rows = []
    for n in sorted(bins):
        d = bins[n]
        bin_rows.append({
            'n_isoforms_bin':        f'{n}+' if n == 10 else str(n),
            'n_genes':               len(d['ufs']) or len(d['ents']),
            'star_mean_uf':          round(mean(d['ufs']), 4)    if d['ufs']    else '',
            'star_median_uf':        round(median(d['ufs']), 4)  if d['ufs']    else '',
            'star_mean_frac_uf0':    round(mean(d['fzeros']), 4) if d['fzeros'] else '',
            'kallisto_mean_entropy': round(mean(d['ents']), 4)   if d['ents']   else '',
            'kallisto_index':        index_label,
            'read_len_bp':           rlen,
            'mode':                  mode,
        })

    bin_fields = ['n_isoforms_bin', 'n_genes', 'star_mean_uf', 'star_median_uf',
                  'star_mean_frac_uf0', 'kallisto_mean_entropy',
                  'kallisto_index', 'read_len_bp', 'mode']
    with open(outdir / 'reads_by_nisoform.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=bin_fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(bin_rows)

    all_ents = [float(r['kallisto_tpm_entropy']) for r in map_rows
                if r['kallisto_tpm_entropy'] != '']
    all_ufs  = [float(r['star_uf_mean']) for r in map_rows if r['star_uf_mean'] != '']
    uf_str   = f'mean_UF={mean(all_ufs):.4f}  ' if all_ufs else '[no STAR UF]  '
    print(f'    [{index_label}]  {uf_str}mean_entropy={mean(all_ents):.4f}')


def _build_mane_kal_lookup(kal_data, gene2ids, enst2sym_conv):
    """
    For a kallisto run against the MANE index (ENST target IDs), rebuild
    the lookup so it is keyed by CDS IDs (NP_) using the id_conversion_table.
    Each gene's CDS entries all map to the same MANE ENST → they receive the
    same per-gene TPM value (no within-gene isoform split in the MANE index).
    """
    # Build enst → gene symbol from conversion table
    # (kal_data keys are ENST; gene2ids keys are gene symbols; need to bridge)
    # Strategy: for each gene, find its ENST in kal_data (strip version if needed)
    kal_by_gene = {}
    for enst, vals in kal_data.items():
        sym = enst2sym_conv.get(enst)
        if sym:
            kal_by_gene[sym] = vals

    # Now re-express as {cds_id: values} so _write_run_outputs can use it
    cds_keyed = {}
    for gene, ids in gene2ids.items():
        if gene in kal_by_gene:
            for cds_id in ids:
                cds_keyed[cds_id] = kal_by_gene[gene]
    return cds_keyed


def analyse_run(label, mode, rlen, read_source, uf_fname, kal_dir,
                gene2ids, meta, enst2sym_conv):
    """
    Process one run condition: write vs_dedup_cds/ and vs_MANE/ outputs.

    enst2sym_conv: {enst_versioned: symbol} from id_conversion_table — used to
    bridge MANE-index (ENST target) kallisto results back to gene symbols.
    """
    check_library(label, rlen, mode)

    run_dir    = HERE / 'runs' / label
    kal_base   = MAPDIR / 'kallisto_results' / kal_dir
    kal_cds_p  = kal_base / 'kallisto_vs_dedup_cds' / 'abundance.tsv'
    kal_mane_p = kal_base / 'kallisto_vs_MANE'      / 'abundance.tsv'

    # STAR selfmap file (CDS-level, only available for dedup_cds PE/SE L100bp)
    star_uf = {}
    if uf_fname:
        uf_path = MAPDIR / 'results' / uf_fname
        if uf_path.exists():
            star_uf = load_star_cds_uf(uf_path)
        else:
            print(f'  [warn] STAR selfmap file not found: {uf_fname}')

    print(f'  {label}')

    # ── vs dedup_cds ──────────────────────────────────────────────────────────
    if kal_cds_p.exists():
        kal_cds = load_kallisto(kal_cds_p)
        _write_run_outputs(run_dir / 'vs_dedup_cds', gene2ids, meta,
                           star_uf, kal_cds, rlen, mode, 'dedup_cds')
    else:
        print(f'    [skip vs_dedup_cds] not found: {kal_cds_p}')

    # ── vs MANE ───────────────────────────────────────────────────────────────
    # The MANE index has one ENST per gene. Re-key by CDS ID so _write_run_outputs
    # can group per-gene entries normally (all CDS isoforms get the gene-level MANE TPM).
    if kal_mane_p.exists():
        kal_mane_raw = load_kallisto(kal_mane_p)
        kal_mane = _build_mane_kal_lookup(kal_mane_raw, gene2ids, enst2sym_conv)
        _write_run_outputs(run_dir / 'vs_MANE', gene2ids, meta,
                           star_uf, kal_mane, rlen, mode, 'MANE')
    else:
        print(f'    [skip vs_MANE] not found: {kal_mane_p}')


def analyse_v8_proper(label, uf_tsv, description, gene2ids, enst2sym_conv):
    """
    v8 proper-pair: MANE RNA reads (150 bp PE, real fragment geometry).
    UF file uses ENST keys — bridge to gene symbols via enst2sym_conv.
    Outputs a gene-level STAR UF table (one row per gene, no isoform split).
    """
    uf_path = MAPDIR / 'results' / uf_tsv
    if not uf_path.exists():
        print(f'  [skip v8_proper] file not found: {uf_tsv}')
        return

    # Load ENST-keyed UF
    enst_uf = {}
    with open(uf_path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            enst_uf[row['transcript_id']] = float(row['uniqueness_factor'])

    # Aggregate to gene level via conversion table
    rows = []
    for gene in sorted(gene2ids):
        ids = gene2ids[gene]
        # Try to find the ENST for this gene in enst2sym_conv (reverse lookup)
        sym2enst_local = {sym: enst for enst, sym in enst2sym_conv.items()
                          if sym == gene}
        enst = next(iter(sym2enst_local.values()), None)
        uf   = enst_uf.get(enst) if enst else None
        rows.append({
            'gene_name':          gene,
            'n_cds_isoforms':     len(ids),
            'enst_id':            enst or '',
            'star_uf':            round(uf, 4) if uf is not None else '',
            'reference':          description,
            'read_source':        'MANE_RNA',
            'read_len_bp':        150,
            'mode':               'PE_proper',
        })

    outdir = HERE / 'runs' / label
    outdir.mkdir(parents=True, exist_ok=True)
    fields = ['gene_name', 'n_cds_isoforms', 'enst_id', 'star_uf',
              'reference', 'read_source', 'read_len_bp', 'mode']
    out = outdir / 'gene_uf.tsv'
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)

    with_uf = [r for r in rows if r['star_uf'] != '']
    ufs = [float(r['star_uf']) for r in with_uf]
    print(f'  {label}  ({description})')
    print(f'    genes with UF: {len(with_uf):,}  mean={mean(ufs):.4f}'
          f'  UF=1.0: {sum(1 for u in ufs if u==1.0):,}'
          f'  UF=0.0: {sum(1 for u in ufs if u==0.0):,}')


def resolve_runs():
    """Check which run conditions have at least one kallisto index available."""
    available = []
    for label, mode, rlen, read_source, uf_fname, kal_dir in _RUNS:
        kal_base = MAPDIR / 'kallisto_results' / kal_dir
        has_cds  = (kal_base / 'kallisto_vs_dedup_cds' / 'abundance.tsv').exists()
        has_mane = (kal_base / 'kallisto_vs_MANE'      / 'abundance.tsv').exists()
        if not (has_cds or has_mane):
            print(f'  [skip] {label} — no kallisto results in {kal_dir}')
            continue
        available.append((label, mode, rlen, read_source, uf_fname, kal_dir))
    return available


# ══════════════════════════════════════════════════════════════════════════════
# E. MAFFT pairwise CDS isoform homology
# ══════════════════════════════════════════════════════════════════════════════

def run_mafft_pairwise(names, sequences):
    """Run MAFFT on sequences and return list of (id1, id2, pct_identity) tuples."""
    fa_in = ''.join(f'>{n}\n{s}\n' for n, s in zip(names, sequences))
    try:
        res = subprocess.run(['mafft', '--quiet', '--auto', '-'],
                             input=fa_in, capture_output=True, text=True, timeout=120)
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return []
    aln, cur = {}, None
    for line in res.stdout.splitlines():
        if line.startswith('>'):
            cur = line[1:].strip(); aln[cur] = []
        elif cur:
            aln[cur].append(line.strip())
    aln  = {k: ''.join(v).upper() for k, v in aln.items()}
    keys = list(aln)
    pairs = []
    for i in range(len(keys)):
        for j in range(i+1, len(keys)):
            a, b = aln[keys[i]], aln[keys[j]]
            if len(a) != len(b): continue
            match   = sum(1 for x, y in zip(a, b) if x == y and x != '-')
            aligned = sum(1 for x, y in zip(a, b) if x != '-' or y != '-')
            pairs.append((keys[i], keys[j],
                          round(100*match/aligned, 2) if aligned else 0))
    return pairs


def run_mafft_cds_isoforms(gene2ids, meta, seqs):
    """Section E: MAFFT on CDS isoforms (within-gene). Cap at MAFFT_MAX_GENES."""
    eligible = sorted(
        [(g, gene2ids[g]) for g in gene2ids
         if len(gene2ids[g]) > 1
         and all(meta[i][3] <= MAFFT_MAX_LEN for i in gene2ids[g])],
        key=lambda x: len(x[1])
    )[:MAFFT_MAX_GENES]

    print(f'  Running MAFFT on {len(eligible)} multi-isoform genes (≤{MAFFT_MAX_LEN} nt) ...')
    hom_rows = []
    for gi, (gene, ids) in enumerate(eligible):
        ns = [i for i in ids if i in seqs]
        ss = [seqs[i] for i in ns]
        if len(ns) < 2: continue
        for id1, id2, pct in run_mafft_pairwise(ns, ss):
            _, p1, _, l1 = meta[id1]
            _, p2, _, l2 = meta[id2]
            hom_rows.append({
                'gene_name':    gene, 'n_isoforms': len(ids),
                'isoform_1':    id1,  'protein_id_1': p1, 'length_nt_1': l1,
                'isoform_2':    id2,  'protein_id_2': p2, 'length_nt_2': l2,
                'pct_identity': pct,  'high_homology': pct >= 95,
            })
        if (gi+1) % 50 == 0:
            print(f'  {gi+1}/{len(eligible)} ...', flush=True)

    fields = ['gene_name', 'n_isoforms', 'isoform_1', 'protein_id_1', 'length_nt_1',
              'isoform_2', 'protein_id_2', 'length_nt_2', 'pct_identity', 'high_homology']
    out = HERE / 'isoform_homology.tsv'
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(hom_rows)
    n_high = sum(1 for r in hom_rows if r['high_homology'])
    print(f'  Saved {out.name}  ({len(hom_rows)} pairs, '
          f'{n_high} [{100*n_high/max(len(hom_rows),1):.1f}%] ≥95% identity)')


# ══════════════════════════════════════════════════════════════════════════════
# F. UTR-fail genes — equal-signal statistical tests
# ══════════════════════════════════════════════════════════════════════════════

def _gene_family_prefix(symbol):
    """
    Extract gene family prefix from HGNC symbol.
    Used to group paralogs for the statistical tests.
    Examples: OR1A1→OR, ZNF123→ZNF, HLA-A→HLA, PCDHA1→PCDHA, HIST1H2AB→HIST1H2A
    """
    # Preserve known multi-segment prefixes
    for prefix in ('HIST1H2A', 'HIST1H2B', 'HIST1H3', 'HIST1H4',
                   'HIST2H2A', 'HIST2H2B', 'HIST3H2A',
                   'PCDHA', 'PCDHB', 'PCDHGA', 'PCDHGB', 'PCDHGC',
                   'DEFB', 'DEFA', 'KRTAP', 'NBPF', 'GOLGA', 'SPDYE', 'PRR20'):
        if symbol.upper().startswith(prefix.upper()):
            return prefix
    # HLA-X → HLA
    if re.match(r'^HLA-', symbol):
        return 'HLA'
    # General: leading uppercase letters up to first digit or hyphen
    m = re.match(r'^([A-Z]+)', symbol)
    return m.group(1) if m else symbol


def _load_utr_fail_genes(kmer_mane_path):
    """
    Load isoform_kmers_mane.tsv and return set of gene names with
    frac_unique < UTR_FAIL_THRESH — these are the UTR-unrescued genes.
    Also returns {gene_name: frac_unique} for all genes.
    """
    fail_set = set()
    frac_by_gene = {}
    if not kmer_mane_path.exists():
        print(f'  [warn] {kmer_mane_path.name} not found — run Section C first.')
        return fail_set, frac_by_gene
    with open(kmer_mane_path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            sym = row['gene_name']
            f   = float(row['frac_unique'])
            frac_by_gene[sym] = f
            if f < UTR_FAIL_THRESH:
                fail_set.add(sym)
    return fail_set, frac_by_gene


def _collect_condition_counts(gene2ids, fail_genes):
    """
    For each UTR-fail gene, collect est_counts across all available conditions.
    Returns {gene: {condition: [est_counts per CDS isoform]}}.
    Only conditions with vs_dedup_cds output are used (CDS isoform split visible).
    """
    data = defaultdict(dict)
    runs_dir = HERE / 'runs'
    for run_dir in sorted(runs_dir.iterdir()):
        cds_map = run_dir / 'vs_dedup_cds' / 'isoform_mappability.tsv'
        if not cds_map.exists():
            continue
        # We need the raw est_counts per CDS isoform, not the per-gene aggregates.
        # Re-load the source kallisto abundance to get per-isoform counts.
        # Infer the kallisto dir from the run label
        label = run_dir.name
        # Map label back to kal_dir
        kal_dir = None
        for entry in _RUNS:
            if entry[0] == label:
                kal_dir = entry[5]
                break
        if kal_dir is None:
            continue
        kal_path = (MAPDIR / 'kallisto_results' / kal_dir
                    / 'kallisto_vs_dedup_cds' / 'abundance.tsv')
        if not kal_path.exists():
            continue
        kal = load_kallisto(kal_path)
        for gene in fail_genes:
            if gene not in gene2ids:
                continue
            ids    = gene2ids[gene]
            counts = [kal[i]['est_counts'] for i in ids if i in kal]
            if len(counts) >= 2:
                data[gene][label] = counts
    return data


def _test_uniform(counts):
    """
    Test if the observed count vector is consistent with uniform distribution.

    Returns a dict with:
      chi2_stat, chi2_p, entropy_norm, cv, n
    chi2 is omitted if scipy unavailable or n < 2.
    """
    n     = len(counts)
    total = sum(counts)
    result = {'n': n, 'total_counts': total}

    if n < 2 or total <= 0:
        result.update({'chi2_stat': '', 'chi2_p': '', 'entropy_norm': 0.0, 'cv': ''})
        return result

    # Normalised Shannon entropy
    props = [c/total for c in counts if c > 0]
    raw_h = -sum(p * log2(p) for p in props)
    result['entropy_norm'] = round(raw_h / log2(n), 4)

    # Coefficient of variation (std/mean across isoforms)
    mean_c = total / n
    result['cv'] = round(stdev(counts) / mean_c, 4) if mean_c > 0 else ''

    # Chi-squared goodness-of-fit against uniform
    if HAS_SCIPY and total > 0:
        expected = [total / n] * n
        try:
            chi2, p = scipy_chisquare(counts, f_exp=expected)
            result['chi2_stat'] = round(chi2, 4)
            result['chi2_p']    = round(p, 6)
        except Exception:
            result['chi2_stat'] = ''
            result['chi2_p']    = ''
    else:
        result['chi2_stat'] = ''
        result['chi2_p']    = ''

    return result


def analyse_utr_fail_stats(gene2ids, frac_by_gene, sym2groups):
    """
    Section F — For each UTR-fail gene, aggregate chi-squared, entropy and CV
    across all 22 conditions (vs_dedup_cds runs only, where the isoform split
    is visible). Genes consistently near-uniform signal are flagged.

    Output: utr_fail_stats.tsv
    Columns: gene, frac_unique_mane, n_cds_isoforms, n_conditions,
             mean_entropy, std_entropy, mean_cv, frac_chi2_p_gt_0.05,
             family_prefix, functional_groups,
             signal_distribution_verdict
    """
    print(f'  Loading UTR-fail gene set (frac_unique < {UTR_FAIL_THRESH}) ...')
    fail_genes = {g for g, f in frac_by_gene.items() if f < UTR_FAIL_THRESH}
    print(f'  {len(fail_genes):,} UTR-fail genes to test')

    cond_counts = _collect_condition_counts(gene2ids, fail_genes)

    rows = []
    for gene in sorted(fail_genes):
        n_iso    = len(gene2ids.get(gene, []))
        fu_mane  = frac_by_gene.get(gene, '')
        prefix   = _gene_family_prefix(gene)
        groups   = ';'.join(sorted(sym2groups.get(gene, []))) or 'none'

        cond_data = cond_counts.get(gene, {})
        if not cond_data:
            rows.append({
                'gene_name':           gene,
                'frac_unique_mane':    fu_mane,
                'n_cds_isoforms':      n_iso,
                'n_conditions':        0,
                'mean_entropy':        '',
                'std_entropy':         '',
                'mean_cv':             '',
                'frac_chi2_nonsig':    '',
                'chi2_note':           'no_data',
                'family_prefix':       prefix,
                'functional_groups':   groups,
                'signal_distribution': 'no_data',
            })
            continue

        stats_per_cond = [_test_uniform(counts)
                          for counts in cond_data.values()]
        entropies = [s['entropy_norm'] for s in stats_per_cond
                     if s['entropy_norm'] != '']
        cvs       = [float(s['cv']) for s in stats_per_cond
                     if s['cv'] not in ('', None)]
        chi2_ps   = [float(s['chi2_p']) for s in stats_per_cond
                     if s['chi2_p'] not in ('', None)]

        mean_ent = round(mean(entropies), 4) if entropies else ''
        std_ent  = round(stdev(entropies), 4) if len(entropies) > 1 else ''
        mean_cv  = round(mean(cvs), 4)        if cvs else ''

        # Fraction of conditions where chi2 p > 0.05 (cannot reject uniform)
        frac_nonsig = round(sum(1 for p in chi2_ps if p > 0.05) / len(chi2_ps), 3) \
                      if chi2_ps else ''
        chi2_note = '' if HAS_SCIPY else 'scipy_missing'

        # Verdict: 'uniform' if consistently high entropy + non-significant chi2
        if mean_ent != '' and mean_ent >= 0.95 and \
                (frac_nonsig == '' or frac_nonsig >= 0.80):
            verdict = 'uniform_signal'      # cannot distinguish isoforms/paralogs
        elif mean_ent != '' and mean_ent >= 0.80:
            verdict = 'mostly_uniform'
        elif mean_ent != '':
            verdict = 'partially_resolved'
        else:
            verdict = 'no_data'

        rows.append({
            'gene_name':           gene,
            'frac_unique_mane':    round(fu_mane, 4) if fu_mane != '' else '',
            'n_cds_isoforms':      n_iso,
            'n_conditions':        len(cond_data),
            'mean_entropy':        mean_ent,
            'std_entropy':         std_ent,
            'mean_cv':             mean_cv,
            'frac_chi2_nonsig':    frac_nonsig,
            'chi2_note':           chi2_note,
            'family_prefix':       prefix,
            'functional_groups':   groups,
            'signal_distribution': verdict,
        })

    fields = ['gene_name', 'frac_unique_mane', 'n_cds_isoforms', 'n_conditions',
              'mean_entropy', 'std_entropy', 'mean_cv', 'frac_chi2_nonsig', 'chi2_note',
              'family_prefix', 'functional_groups', 'signal_distribution']
    out = HERE / 'utr_fail_stats.tsv'
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)

    n_uniform = sum(1 for r in rows if r['signal_distribution'] == 'uniform_signal')
    n_mostly  = sum(1 for r in rows if r['signal_distribution'] == 'mostly_uniform')
    print(f'  Saved {out.name}  ({len(rows):,} genes)')
    print(f'  uniform_signal: {n_uniform:,}  mostly_uniform: {n_mostly:,}  '
          f'partially_resolved: {len(rows)-n_uniform-n_mostly:,}')
    print(f'\n  NOTE: Section F ANOVA per functional group is a separate follow-up.')
    print(f'  Columns mean_entropy, mean_cv, and frac_chi2_nonsig are the inputs.')
    print(f'  Run: R> aov(mean_entropy ~ functional_groups, data=utr_fail_stats)')
    return rows


# ══════════════════════════════════════════════════════════════════════════════
# G. Paralog family MAFFT (Section G)
# ══════════════════════════════════════════════════════════════════════════════

def _build_paralog_families(fail_genes, gene2ids, gencode_symbol2ensts):
    """
    Group UTR-fail genes into paralog families by gene name prefix.
    Filter: 2 ≤ family size ≤ PARALOG_MAX.
    Returns: {prefix: [gene_symbol, ...]}
    """
    prefix_map = defaultdict(list)
    for gene in fail_genes:
        if gene in gene2ids or gene in gencode_symbol2ensts:
            prefix_map[_gene_family_prefix(gene)].append(gene)

    families = {p: sorted(members)
                for p, members in prefix_map.items()
                if 2 <= len(members) <= PARALOG_MAX}
    return families


def _verify_ids(families, gene2ids, symbol2enst_mane, gencode_symbol2ensts,
                enst2sym_conv, np2sym):
    """
    Cross-check IDs for all paralog family members across four sources:
      1. CDS FASTA (gene2ids — NP_ IDs)
      2. MANE RNA FASTA (symbol2enst_mane — ENST IDs)
      3. GENCODE pc_transcripts (gencode_symbol2ensts — ENST IDs)
      4. id_conversion_table (np2sym — NP_ → symbol)

    Writes id_verification.tsv with one row per gene × source.
    """
    all_genes = sorted({g for members in families.values() for g in members})
    rows = []
    for gene in all_genes:
        rows.append({
            'gene_name':              gene,
            'family_prefix':          _gene_family_prefix(gene),
            'in_CDS_FASTA':           int(gene in gene2ids),
            'n_CDS_isoforms':         len(gene2ids.get(gene, [])),
            'MANE_ENST':              symbol2enst_mane.get(gene, ''),
            'in_MANE_FASTA':          int(gene in symbol2enst_mane),
            'n_GENCODE_isoforms':     len(gencode_symbol2ensts.get(gene, [])),
            'in_GENCODE_FA':          int(gene in gencode_symbol2ensts),
            'conv_table_symbol_ok':   int(any(
                                          np2sym.get(np_) == gene
                                          for np_ in gene2ids.get(gene, [])
                                          if not np_.startswith('lcl'))),
            'symbol_consistent':      int(
                                          (gene in gene2ids or gene in symbol2enst_mane)
                                          and (not gene2ids.get(gene) or
                                               any(np2sym.get(np_) == gene
                                                   for np_ in gene2ids.get(gene, [])
                                                   if not np_.startswith('lcl'))
                                               or gene in symbol2enst_mane)),
        })

    out = HERE / 'id_verification.tsv'
    fields = ['gene_name', 'family_prefix', 'in_CDS_FASTA', 'n_CDS_isoforms',
              'MANE_ENST', 'in_MANE_FASTA', 'n_GENCODE_isoforms', 'in_GENCODE_FA',
              'conv_table_symbol_ok', 'symbol_consistent']
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)

    inconsistent = [r for r in rows if not r['symbol_consistent']]
    print(f'  id_verification.tsv: {len(rows)} genes, '
          f'{len(inconsistent)} with symbol inconsistencies')
    if inconsistent:
        for r in inconsistent[:10]:
            print(f'    [!] {r["gene_name"]}: CDS={r["in_CDS_FASTA"]} '
                  f'MANE={r["in_MANE_FASTA"]} GENCODE={r["in_GENCODE_FA"]} '
                  f'conv_ok={r["conv_table_symbol_ok"]}')
    return out


def run_paralog_mafft(families, symbol2enst_mane, enst2seq_mane,
                      gencode_symbol2ensts, enst2seq_gencode, gene2ids):
    """
    Section G: MAFFT on paralog families (cross-gene, not within-gene isoforms).

    Option A: one MANE canonical sequence per family member (full UTR-CDS-UTR).
              Source: MANE RNA FASTA.
    Option B: all GENCODE isoforms per family member if n_isoforms ≤ PARALOG_ISOFORM_MAX.
              Source: GENCODE v49 pc_transcripts.fa.
    """
    rows_a, rows_b = [], []
    total = len(families)

    for fi, (prefix, members) in enumerate(sorted(families.items())):
        print(f'  [{fi+1}/{total}] {prefix}  ({len(members)} members)', flush=True)

        # ── Option A: MANE canonical ─────────────────────────────────────────
        names_a, seqs_a = [], []
        for gene in members:
            enst = symbol2enst_mane.get(gene)
            if enst and enst in enst2seq_mane:
                names_a.append(gene)
                seqs_a.append(enst2seq_mane[enst])

        if len(names_a) >= 2:
            for g1, g2, pct in run_mafft_pairwise(names_a, seqs_a):
                rows_a.append({
                    'family_prefix': prefix,
                    'n_family':      len(members),
                    'gene_1':        g1,
                    'gene_2':        g2,
                    'pct_identity':  pct,
                    'high_homology': pct >= 85,
                    'option':        'A_MANE_canonical',
                })

        # ── Option B: all GENCODE isoforms per member ─────────────────────────
        # Only attempt if every member has ≤ PARALOG_ISOFORM_MAX isoforms.
        isoform_counts = {g: len(gencode_symbol2ensts.get(g, []))
                          for g in members}
        max_iso = max(isoform_counts.values()) if isoform_counts else 0

        if max_iso <= PARALOG_ISOFORM_MAX and max_iso > 0:
            names_b, seqs_b = [], []
            for gene in members:
                ensts = gencode_symbol2ensts.get(gene, [])
                for i, enst in enumerate(ensts):
                    if enst in enst2seq_gencode:
                        names_b.append(f'{gene}_iso{i+1}')
                        seqs_b.append(enst2seq_gencode[enst])

            if len(names_b) >= 2:
                for id1, id2, pct in run_mafft_pairwise(names_b, seqs_b):
                    g1 = id1.rsplit('_iso', 1)[0]
                    g2 = id2.rsplit('_iso', 1)[0]
                    rows_b.append({
                        'family_prefix':    prefix,
                        'n_family':         len(members),
                        'gene_1':           g1,
                        'gene_2':           g2,
                        'isoform_label_1':  id1,
                        'isoform_label_2':  id2,
                        'pct_identity':     pct,
                        'high_homology':    pct >= 85,
                        'option':           'B_GENCODE_all_isoforms',
                    })
        else:
            print(f'    [skip B] {prefix}: max isoforms/member={max_iso} '
                  f'(limit {PARALOG_ISOFORM_MAX})')

    # Write Option A
    out_a = HERE / 'paralog_homology_mane.tsv'
    if rows_a:
        fields_a = ['family_prefix', 'n_family', 'gene_1', 'gene_2',
                    'pct_identity', 'high_homology', 'option']
        with open(out_a, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=fields_a, delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(rows_a)
        n_high_a = sum(1 for r in rows_a if r['high_homology'])
        print(f'  Saved {out_a.name}  ({len(rows_a)} pairs, '
              f'{n_high_a} [{100*n_high_a/max(len(rows_a),1):.1f}%] ≥85% identity)')
    else:
        print(f'  [warn] no Option A results — check MANE FASTA coverage')

    # Write Option B
    out_b = HERE / 'paralog_homology_isoforms.tsv'
    if rows_b:
        fields_b = ['family_prefix', 'n_family', 'gene_1', 'gene_2',
                    'isoform_label_1', 'isoform_label_2',
                    'pct_identity', 'high_homology', 'option']
        with open(out_b, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=fields_b, delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(rows_b)
        n_high_b = sum(1 for r in rows_b if r['high_homology'])
        print(f'  Saved {out_b.name}  ({len(rows_b)} pairs, '
              f'{n_high_b} [{100*n_high_b/max(len(rows_b),1):.1f}%] ≥85% identity)')
    else:
        print(f'  [info] no Option B results (all families exceed {PARALOG_ISOFORM_MAX} '
              f'isoforms/member or no GENCODE coverage)')

    return rows_a, rows_b


# ══════════════════════════════════════════════════════════════════════════════
# H. Actionable gene lists
# ══════════════════════════════════════════════════════════════════════════════

def build_actionable_lists(utr_fail_rows, rows_a, rows_b, frac_by_gene, gene2ids):
    """
    Section H — Two actionable output lists:

    gene_level_recommendation.tsv
        Genes where isoform redundancy makes the CDS isoform signal indistinguishable:
        high kallisto entropy + uniform chi-squared signal + low UTR rescue.
        Recommendation: collapse all isoforms to gene-level TPM sum.

    cluster_level_recommendation.tsv
        Genes in paralog families where pairwise identity ≥ 85% (MAFFT Option A)
        AND UTR rescue is insufficient (frac_unique < UTR_FAIL_THRESH).
        Recommendation: aggregate all family members into a paralog-cluster sum;
        do not interpret individual gene TPMs.
    """

    # ── Gene-level: isoform redundancy ───────────────────────────────────────
    gene_rows = []
    for r in utr_fail_rows:
        if r['n_cds_isoforms'] < 2:
            continue
        if r['signal_distribution'] not in ('uniform_signal', 'mostly_uniform'):
            continue
        gene_rows.append({
            'gene_name':               r['gene_name'],
            'n_cds_isoforms':          r['n_cds_isoforms'],
            'frac_unique_mane':        r['frac_unique_mane'],
            'mean_entropy':            r['mean_entropy'],
            'mean_cv':                 r['mean_cv'],
            'frac_chi2_nonsig':        r['frac_chi2_nonsig'],
            'signal_distribution':     r['signal_distribution'],
            'functional_groups':       r['functional_groups'],
            'recommendation':          'collapse_to_gene_level',
            'rationale':               ('Multiple CDS isoforms with high entropy '
                                        'signal — kallisto cannot distinguish isoforms. '
                                        'UTR k-mers provide insufficient rescue '
                                        f'(frac_unique={r["frac_unique_mane"]}).'),
        })

    gene_fields = ['gene_name', 'n_cds_isoforms', 'frac_unique_mane',
                   'mean_entropy', 'mean_cv', 'frac_chi2_nonsig',
                   'signal_distribution', 'functional_groups',
                   'recommendation', 'rationale']
    out_g = HERE / 'gene_level_recommendation.tsv'
    with open(out_g, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=gene_fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(gene_rows)
    print(f'  Saved {out_g.name}  ({len(gene_rows):,} genes → gene-level)')

    # ── Cluster-level: paralog ambiguity ─────────────────────────────────────
    # Build {gene: max pairwise identity with any family member} from Option A
    gene_max_ident = defaultdict(float)
    gene_cluster   = {}
    for r in rows_a:
        if r['pct_identity'] > gene_max_ident[r['gene_1']]:
            gene_max_ident[r['gene_1']] = r['pct_identity']
            gene_cluster[r['gene_1']]   = r['family_prefix']
        if r['pct_identity'] > gene_max_ident[r['gene_2']]:
            gene_max_ident[r['gene_2']] = r['pct_identity']
            gene_cluster[r['gene_2']]   = r['family_prefix']

    cluster_rows = []
    for gene, max_ident in sorted(gene_max_ident.items()):
        fu = frac_by_gene.get(gene, 1.0)
        if max_ident < 85.0:
            continue
        if fu >= UTR_FAIL_THRESH:
            continue  # UTR rescue works for this gene
        cluster_rows.append({
            'gene_name':           gene,
            'family_prefix':       gene_cluster.get(gene, ''),
            'max_paralog_identity': max_ident,
            'frac_unique_mane':    round(fu, 4),
            'n_cds_isoforms':      len(gene2ids.get(gene, [])),
            'recommendation':      'aggregate_at_cluster_level',
            'rationale':           (f'Paralog identity ≥ {max_ident:.1f}% with at least '
                                    'one family member AND UTR k-mers fail to rescue '
                                    f'(frac_unique={round(fu, 3)}). '
                                    'Individual gene TPMs are unreliable; use cluster sum.'),
        })

    cluster_fields = ['gene_name', 'family_prefix', 'max_paralog_identity',
                      'frac_unique_mane', 'n_cds_isoforms',
                      'recommendation', 'rationale']
    out_c = HERE / 'cluster_level_recommendation.tsv'
    with open(out_c, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cluster_fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(cluster_rows)
    print(f'  Saved {out_c.name}  ({len(cluster_rows):,} genes → cluster-level)')

    # Overlap summary
    gene_set    = {r['gene_name'] for r in gene_rows}
    cluster_set = {r['gene_name'] for r in cluster_rows}
    both        = gene_set & cluster_set
    print(f'\n  Overlap (isoform AND paralog ambiguity): {len(both)} genes')
    if both:
        print(f'  Examples: {", ".join(sorted(both)[:10])}')


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    print('=' * 72)
    print('01_isoform_redundancy  [v2 — full 22-condition + UTR-fail + paralog MAFFT]')
    print('=' * 72)

    # ── A. Parse CDS FASTA ────────────────────────────────────────────────────
    print('\n[A] Parsing dedup CDS FASTA ...')
    gene2ids, meta, seqs = parse_cds_fasta(CDS_FA)
    print(f'  {len(seqs):,} sequences  |  {len(gene2ids):,} genes')
    n1 = sum(1 for ids in gene2ids.values() if len(ids) == 1)
    nm = len(gene2ids) - n1
    print(f'  single-isoform: {n1:,} ({100*n1/len(gene2ids):.1f}%)  '
          f'multi-isoform: {nm:,} ({100*nm/len(gene2ids):.1f}%)')
    iso_rows = make_isoform_table(gene2ids, meta)
    write_isoform_table(iso_rows)

    # ── B. Functional groups + pseudogenes ────────────────────────────────────
    print('\n[B] Functional groups and pseudogene counts ...')
    np2sym, enst2sym_conv = load_id_conversion(CONV_TSV)
    sym2groups = load_groups(enst2sym_conv)
    pseudo_counts = count_pseudogenes_gencode(set(gene2ids.keys()))
    group_rows = write_isoform_groups(
        iso_rows, pseudo_counts, sym2groups, np2sym, meta, gene2ids)

    # ── C. 31-mer unique fraction + unique stretch ────────────────────────────
    print(f'\n[C] {KMER_K}-mer unique fraction ...')
    print('  CDS within-gene (multi-isoform genes only) ...')
    cds_kmer_rows = kmer_unique_fraction_within_gene(gene2ids, seqs, KMER_K)
    write_kmer_table(cds_kmer_rows, HERE / 'isoform_kmers_cds.tsv')

    print('  MANE RNA cross-gene (UTR + junction k-mers) ...')
    symbol2enst_mane, enst2seq_mane = {}, {}
    if MANE_FA.exists():
        symbol2enst_mane, enst2seq_mane = parse_mane_fasta(MANE_FA)
        mane_kmer_rows = kmer_unique_fraction_mane(symbol2enst_mane, enst2seq_mane, KMER_K)
        write_kmer_table(mane_kmer_rows, HERE / 'isoform_kmers_mane.tsv')
        # Build gene_name → frac_unique map for downstream sections
        frac_by_gene = {r['gene_name']: r['frac_unique'] for r in mane_kmer_rows}
    else:
        print(f'  [skip] MANE RNA FASTA not found: {MANE_FA}')
        mane_kmer_rows = []
        frac_by_gene = {}

    # ── D. Per-run variants (22 conditions × 2 indices + v8 proper) ───────────
    print('\n[D] Run variants ...')
    runs = resolve_runs()
    print(f'  {len(runs)} available conditions')
    for label, mode, rlen, read_source, uf_fname, kal_dir in runs:
        analyse_run(label, mode, rlen, read_source, uf_fname, kal_dir,
                    gene2ids, meta, enst2sym_conv)

    print('\n  v8 proper-pair runs (ENST-keyed, 150 bp PE proper geometry):')
    for label, uf_tsv, description in _V8_PROPER:
        analyse_v8_proper(label, uf_tsv, description, gene2ids, enst2sym_conv)

    # ── E. MAFFT CDS isoform homology ─────────────────────────────────────────
    print('\n[E] MAFFT pairwise CDS isoform homology ...')
    run_mafft_cds_isoforms(gene2ids, meta, seqs)

    # ── F. UTR-fail statistical tests ─────────────────────────────────────────
    print(f'\n[F] UTR-fail equal-signal tests (χ², entropy, CV) ...')
    utr_fail_rows = []
    if frac_by_gene:
        utr_fail_rows = analyse_utr_fail_stats(gene2ids, frac_by_gene, sym2groups)
    else:
        print('  [skip] isoform_kmers_mane.tsv not available')

    # ── G. Paralog MAFFT (families 2–99 members) ─────────────────────────────
    print(f'\n[G] Paralog family MAFFT (≤{PARALOG_MAX} members) ...')
    rows_a, rows_b = [], []
    if frac_by_gene and GENCODE_FA.exists():
        print('  Parsing GENCODE v49 protein-coding transcripts ...')
        gencode_sym2ensts, enst2seq_gencode, _ = parse_gencode_pc_fasta(GENCODE_FA)
        print(f'  {len(gencode_sym2ensts):,} genes, '
              f'{len(enst2seq_gencode):,} isoform sequences')

        fail_genes = {g for g, f in frac_by_gene.items() if f < UTR_FAIL_THRESH}
        families   = _build_paralog_families(
            fail_genes, gene2ids, gencode_sym2ensts)
        print(f'  {len(families)} paralog families identified '
              f'(prefix-based, 2–{PARALOG_MAX} members)')

        # ID cross-verification before running MAFFT
        print('  Verifying IDs across CDS/MANE/GENCODE/conversion table ...')
        _verify_ids(families, gene2ids, symbol2enst_mane,
                    gencode_sym2ensts, enst2sym_conv, np2sym)

        rows_a, rows_b = run_paralog_mafft(
            families, symbol2enst_mane, enst2seq_mane,
            gencode_sym2ensts, enst2seq_gencode, gene2ids)
    elif not GENCODE_FA.exists():
        print(f'  [skip] GENCODE pc_transcripts.fa not found: {GENCODE_FA}')
    else:
        print('  [skip] no frac_by_gene data (Section C must run first)')

    # ── H. Actionable gene lists ───────────────────────────────────────────────
    print('\n[H] Building actionable gene lists ...')
    if utr_fail_rows or rows_a:
        build_actionable_lists(utr_fail_rows, rows_a, rows_b, frac_by_gene, gene2ids)
    else:
        print('  [skip] Sections F and G produced no data')

    print('\n' + '=' * 72)
    print('Done.')


if __name__ == '__main__':
    main()
