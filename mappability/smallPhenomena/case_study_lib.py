"""
case_study_lib.py — Shared functions for 08A / 08B gene-level case studies.

Provides:
  build_enst_gene_map(fasta_path)
      Parse MANE RNA FASTA headers → {ENST_no_version: gene_symbol}

  load_gene_annotation(annot_path)
      Load gene_overlap_annotation.tsv → {gene_name: row_dict}
      (first protein-coding row wins if multiple entries for same name)

  load_all_star_runs(sc3utr_dir)
      Scan sc3UTRs/results/ for all sc3prime_*.tsv files.
      Returns {label: {ENST_no_version: row_dict}}
      Includes read_len and mode decoded from the label.

  load_all_kallisto_runs(kal_dir)
      Scan kallisto_results/ for all mappability_MANE_* subdirectories.
      Returns {run_name: {ENST_no_version: {'tpm': float, 'est_counts': float}}}

  classify_gene(annot_row, gene_symbol)
      Returns list of ambiguity category strings for the gene.

  STAR_RUN_META   dict: label → (read_len_bp, mode, read_source, reference)
  KAL_RUN_META    dict: run_name → (read_len_bp, mode, read_source, reference)
"""

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Known run metadata so we can annotate every result row
# ---------------------------------------------------------------------------

# sc3prime label → (read_len_bp, mode, read_source, reference)
STAR_RUN_META = {
    # MANE RNA
    'v5a_MANE_RNA_SE_100':     (100, 'SE', 'RNA', 'MANE'),
    'v5c_MANE_RNA_SE_75':      ( 75, 'SE', 'RNA', 'MANE'),
    'v5e_MANE_RNA_SE_150':     (150, 'SE', 'RNA', 'MANE'),
    'v5d_MANE_RNA_SE_200':     (200, 'SE', 'RNA', 'MANE'),
    'v5f_MANE_RNA_PE_75':      ( 75, 'PE', 'RNA', 'MANE'),
    'v5g_MANE_RNA_PE_100':     (100, 'PE', 'RNA', 'MANE'),
    'v5h_MANE_RNA_PE_150':     (150, 'PE', 'RNA', 'MANE'),
    'v5i_MANE_RNA_PE_200':     (200, 'PE', 'RNA', 'MANE'),
    # MANE exon-flank reads (flanked simulation: wider context)
    'v5b_MANE_exflank_PE_100': (100, 'PE', 'exflank', 'MANE'),
    'v5j_MANE_exflank_SE_75':  ( 75, 'SE', 'exflank', 'MANE'),
    'v5k_MANE_exflank_SE_100': (100, 'SE', 'exflank', 'MANE'),
    'v5l_MANE_exflank_SE_150': (150, 'SE', 'exflank', 'MANE'),
    'v5m_MANE_exflank_PE_75':  ( 75, 'PE', 'exflank', 'MANE'),
    'v5n_MANE_exflank_PE_150': (150, 'PE', 'exflank', 'MANE'),
    # Proper paired-end
    'v8c_MANE_RNA_PE_150_proper': (150, 'PE_proper', 'RNA', 'MANE'),
    # Full-genome RNA (for comparison)
    'v1c_genome_RNA_SE_100':   (100, 'SE', 'RNA', 'genome'),
    'v1d_genome_RNA_PE_100':   (100, 'PE', 'RNA', 'genome'),
    'v3c_genome_RNA_SE_75':    ( 75, 'SE', 'RNA', 'genome'),
    'v3d_genome_RNA_PE_75':    ( 75, 'PE', 'RNA', 'genome'),
    'v4c_genome_RNA_SE_150':   (150, 'SE', 'RNA', 'genome'),
    'v4d_genome_RNA_PE_150':   (150, 'PE', 'RNA', 'genome'),
    'v6a_genome_RNA_SE_200':   (200, 'SE', 'RNA', 'genome'),
    'v6b_genome_RNA_PE_200':   (200, 'PE', 'RNA', 'genome'),
    'v8a_genome_RNA_PE_150_proper': (150, 'PE_proper', 'RNA', 'genome'),
    # dedup CDS reference
    'v2c_dedupcds_RNA_SE_100': (100, 'SE', 'RNA', 'dedup_cds'),
    'v2e_dedupcds_RNA_PE_100': (100, 'PE', 'RNA', 'dedup_cds'),
    'v2j_dedupcds_RNA_SE_75':  ( 75, 'SE', 'RNA', 'dedup_cds'),
    'v2k_dedupcds_RNA_PE_75':  ( 75, 'PE', 'RNA', 'dedup_cds'),
    'v2l_dedupcds_RNA_SE_150': (150, 'SE', 'RNA', 'dedup_cds'),
    'v2m_dedupcds_RNA_PE_150': (150, 'PE', 'RNA', 'dedup_cds'),
    'v2n_dedupcds_RNA_SE_200': (200, 'SE', 'RNA', 'dedup_cds'),
    'v2o_dedupcds_RNA_PE_200': (200, 'PE', 'RNA', 'dedup_cds'),
    'v8b_dedupcds_RNA_PE_150_proper': (150, 'PE_proper', 'RNA', 'dedup_cds'),
}

# kallisto run dir suffix → (read_len_bp, mode, read_source, reference_used)
KAL_RUN_META = {
    'mappability_MANE_RNA_PE':         (100, 'PE',  'RNA',     'MANE'),
    'mappability_MANE_RNA_PE_L75bp':   ( 75, 'PE',  'RNA',     'MANE'),
    'mappability_MANE_RNA_PE_L150bp':  (150, 'PE',  'RNA',     'MANE'),
    'mappability_MANE_RNA_PE_L200bp':  (200, 'PE',  'RNA',     'MANE'),
    'mappability_MANE_RNA_SE':         (100, 'SE',  'RNA',     'MANE'),
    'mappability_MANE_RNA_SE_L75bp':   ( 75, 'SE',  'RNA',     'MANE'),
    'mappability_MANE_RNA_SE_L150bp':  (150, 'SE',  'RNA',     'MANE'),
    'mappability_MANE_RNA_SE_L200bp':  (200, 'SE',  'RNA',     'MANE'),
    'mappability_MANE_exonflank_PE':   (100, 'PE',  'exflank', 'MANE'),
    'mappability_MANE_exonflank_PE_L75bp':  ( 75, 'PE', 'exflank', 'MANE'),
    'mappability_MANE_exonflank_PE_L150bp': (150, 'PE', 'exflank', 'MANE'),
    'mappability_MANE_exonflank_SE':   (100, 'SE',  'exflank', 'MANE'),
    'mappability_MANE_exonflank_SE_L75bp':  ( 75, 'SE', 'exflank', 'MANE'),
    'mappability_MANE_exonflank_SE_L150bp': (150, 'SE', 'exflank', 'MANE'),
}

# Ambiguity category labels (priority-ordered)
AMBIGUITY_CATEGORIES = [
    'isoform_ambiguity',       # multiple CDS isoforms in the reference
    'paralog_ambiguity',       # multiple coding copies with high homology
    'local_duplication',       # tandem/proximal duplicate copies
    'pseudogene_ambiguity',    # co-located with processed/unprocessed pseudogene
    'genomic_overlap',         # overlapping annotation on same locus
    'nested_gene',             # one gene nested inside another
    'allelic_ambiguity',       # high allelic diversity (HLA, highly polymorphic)
    'multi_copy_identical',    # multiple nearly-identical loci (RNU, histones)
    'repeat_ambiguity',        # inside or flanked by repeat elements (LINEs, HERVs)
    'readthrough',             # GENCODE readthrough_gene
    'retrogene',               # retrogene insertion
    'unique',                  # no known ambiguity source
]

# Genes with known high allelic diversity (not derivable from annotation alone)
HLA_CLASS_I_GENES  = {'HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G',
                       'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-T', 'HLA-V'}
HLA_CLASS_II_GENES = {'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5',
                       'HLA-DRB6', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DQA2', 'HLA-DQB2',
                       'HLA-DPA1', 'HLA-DPB1', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB'}
HLA_GENES = HLA_CLASS_I_GENES | HLA_CLASS_II_GENES


def build_enst_gene_map(fasta_path):
    """
    Parse MANE RNA FASTA headers.
    Format: >ENST... cdna ... gene_symbol:SYMBOL ...
    Returns {ENST_no_version: gene_symbol}
    """
    mapping = {}
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            enst = line[1:].split()[0].split('.')[0]
            m = re.search(r'gene_symbol:(\S+)', line)
            mapping[enst] = m.group(1) if m else ''
    return mapping


def load_gene_annotation(annot_path):
    """
    Returns two dicts:
      by_name  {gene_name: annot_row_dict}  — protein-coding first, else first
      by_id    {gene_id: annot_row_dict}
    """
    by_name = {}
    by_id   = {}
    with open(annot_path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            gid  = row['gene_id']
            name = row['gene_name']
            by_id[gid] = row
            # protein-coding entry wins over any previously stored entry
            if name not in by_name or row['gene_type'] == 'protein_coding':
                by_name[name] = row
    return by_name, by_id


def load_all_star_runs(sc3utr_dir, mane_only=True):
    """
    Scan sc3UTRs/results/ for sc3prime_*.tsv files.
    Returns {label: {ENST_no_version: row_dict}}
    If mane_only=True, only loads files whose label is in STAR_RUN_META
    AND whose reference is 'MANE'.
    """
    results = {}
    sc3dir  = Path(sc3utr_dir)
    for tsv in sorted(sc3dir.glob('sc3prime_*.tsv')):
        # Extract label from first data row
        with open(tsv) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            first  = next(iter(reader), None)
            if first is None:
                continue
            label = first.get('label', '')
        if label not in STAR_RUN_META:
            continue
        meta = STAR_RUN_META[label]
        if mane_only and meta[3] != 'MANE':
            continue
        # Reload full file
        run_data = {}
        with open(tsv) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                enst = row['transcript_id'].split('.')[0]
                run_data[enst] = row
        results[label] = run_data
    return results


def load_all_kallisto_runs(kal_dir):
    """
    Scan kallisto_results/ for mappability_MANE_* directories.
    Returns {run_name: {ENST_no_version: {'tpm': float, 'est_counts': float}}}
    NOTE: est_counts are NOT library-size normalised — do not compare across runs.
          TPM is per-million within each run — cross-sample fraction comparison OK.
    """
    results = {}
    kal_path = Path(kal_dir)
    for run_dir in sorted(kal_path.iterdir()):
        name = run_dir.name
        if name not in KAL_RUN_META:
            continue
        # Try kallisto_vs_MANE first, then kallisto_vs_dedup_cds
        for sub in ('kallisto_vs_MANE', 'kallisto_vs_dedup_cds'):
            abund = run_dir / sub / 'abundance.tsv'
            if abund.exists():
                run_data = {}
                with open(abund) as fh:
                    for row in csv.DictReader(fh, delimiter='\t'):
                        enst = row['target_id'].split('.')[0]
                        run_data[enst] = {
                            'tpm':        float(row['tpm']),
                            'est_counts': float(row['est_counts']),
                        }
                results[name] = run_data
                break
    return results


def classify_gene(annot_row, gene_symbol):
    """
    Returns sorted list of ambiguity category strings for the gene.
    Uses the annot_row from gene_overlap_annotation.tsv.
    """
    categories = set()

    if annot_row is None:
        return ['not_in_annotation']

    phen = annot_row.get('phenomenon', '')
    gt   = annot_row.get('gene_type', '')

    # Primary phenomenon mapping
    phen_map = {
        'genomic_overlap_antisense': 'genomic_overlap',
        'genomic_overlap_sense':     'paralog_ambiguity',
        'nested':                    'nested_gene',
        'local_tandem_dup':          'local_duplication',
        'distant_paralog':           'paralog_ambiguity',
        'pseudogene_shadow':         'pseudogene_ambiguity',
        'readthrough':               'readthrough',
        'retrogene':                 'retrogene',
        'unique':                    'unique',
    }
    if phen in phen_map:
        categories.add(phen_map[phen])

    # Sequence cluster with multiple PC members → paralog/multi-copy
    try:
        n_pc_cluster = int(annot_row.get('seq_cluster_n_pc', 0))
    except ValueError:
        n_pc_cluster = 0
    if n_pc_cluster > 2:
        categories.add('paralog_ambiguity')
    elif n_pc_cluster == 2:
        # Check if it's truly a known multi-copy (same chromosome, local)
        try:
            span = float(annot_row.get('seq_cluster_max_span_mb', 999))
        except ValueError:
            span = 999
        multi_chrom = annot_row.get('seq_cluster_multi_chrom', 'True').lower()
        if span <= 5 and multi_chrom in ('false', '0', 'no'):
            categories.add('local_duplication')
        else:
            categories.add('paralog_ambiguity')

    # Pseudogene in overlap group
    biotypes = annot_row.get('overlap_group_biotypes', '')
    if any(b in biotypes for b in ['pseudogene', 'processed_pseudogene',
                                    'unprocessed_pseudogene']):
        categories.add('pseudogene_ambiguity')

    # Isoform ambiguity: if gene has multiple CDS isoforms
    # (detected from seq_cluster_size when mixed=False, PC only)
    try:
        cluster_sz = int(annot_row.get('seq_cluster_size', 0))
        cluster_npseudo = int(annot_row.get('seq_cluster_n_pseudo', 0))
    except ValueError:
        cluster_sz = cluster_npseudo = 0
    if cluster_sz > 1 and n_pc_cluster <= 1 and cluster_npseudo == 0:
        categories.add('isoform_ambiguity')

    # Allelic ambiguity: HLA genes and other known highly polymorphic loci
    if gene_symbol in HLA_GENES or gene_symbol.startswith('HLA-'):
        categories.add('allelic_ambiguity')

    # Multi-copy identical: small ncRNA genes with many copies (RNU, SNRP…)
    if gt in ('snRNA', 'snoRNA') or re.match(r'RNU\d|SNRP', gene_symbol):
        categories.add('multi_copy_identical')

    if not categories or categories == {'unique'}:
        categories = {'unique'}

    return sorted(categories)


def extract_gene_rows(target_genes, enst_gene_map, gene_annot_by_name,
                      star_runs, kal_runs):
    """
    For each gene in target_genes, return:
      (1) classification dict
      (2) list of per-run STAR rows
      (3) list of per-run kallisto rows

    enst_gene_map   : {ENST: gene_symbol}  — built from MANE FASTA
    gene_annot_by_name : {gene_name: annot_row}
    star_runs       : {label: {ENST: row_dict}}
    kal_runs        : {run_name: {ENST: {tpm, est_counts}}}
    """
    # Reverse map: gene_symbol → list of ENST IDs
    gene_to_enst = defaultdict(list)
    for enst, sym in enst_gene_map.items():
        if sym:
            gene_to_enst[sym].append(enst)

    # Infer ENST for each target gene
    classifications = {}
    star_detail     = []
    kal_detail      = []
    missing_genes   = []

    for gsym in target_genes:
        enst_ids = gene_to_enst.get(gsym, [])
        arow     = gene_annot_by_name.get(gsym)
        cats     = classify_gene(arow, gsym)
        in_mane  = len(enst_ids) > 0

        classifications[gsym] = {
            'gene_symbol':       gsym,
            'in_mane':           in_mane,
            'n_enst_ids':        len(enst_ids),
            'enst_ids':          ';'.join(enst_ids),
            'chrom':             arow.get('chrom', '') if arow else 'N/A',
            'strand':            arow.get('strand', '') if arow else 'N/A',
            'gene_type':         arow.get('gene_type', '') if arow else 'N/A',
            'phenomenon':        arow.get('phenomenon', '') if arow else 'N/A',
            'overlap_group_size': arow.get('overlap_group_size', '') if arow else '',
            'overlap_group_biotypes': arow.get('overlap_group_biotypes', '') if arow else '',
            'overlap_group_genes': arow.get('overlap_group_genes', '') if arow else '',
            'seq_cluster_n_pc':  arow.get('seq_cluster_n_pc', '') if arow else '',
            'seq_cluster_max_span_mb': arow.get('seq_cluster_max_span_mb', '') if arow else '',
            'ambiguity_types':   ';'.join(cats),
            'is_HLA':            gsym in HLA_GENES or gsym.startswith('HLA-'),
        }
        if not in_mane:
            missing_genes.append(gsym)

        # STAR rows — one per run × ENST
        for label, run_data in sorted(star_runs.items()):
            meta = STAR_RUN_META.get(label, (None, None, None, None))
            for enst in enst_ids:
                row = run_data.get(enst, {})
                star_detail.append({
                    'gene_symbol':   gsym,
                    'enst_id':       enst,
                    'run_label':     label,
                    'read_len_bp':   meta[0],
                    'mode':          meta[1],
                    'read_source':   meta[2],
                    'reference':     meta[3],
                    'transcript_len_nt': row.get('transcript_len', ''),
                    'uf_overall':    row.get('uf_overall', ''),
                    'uf_3p250':      row.get('uf_3p250', ''),
                    'delta_3p250':   row.get('delta_3p250', ''),
                    'uf_3p400':      row.get('uf_3p400', ''),
                    'delta_3p400':   row.get('delta_3p400', ''),
                    'uf_3p500':      row.get('uf_3p500', ''),
                    'delta_3p500':   row.get('delta_3p500', ''),
                    'n_positions':   row.get('n_positions', ''),
                    'n_unique_3p250': row.get('n_unique_3p250', ''),
                    'ambiguity_types': ';'.join(cats),
                })

        # Kallisto rows — one per run × ENST
        for run_name, run_data in sorted(kal_runs.items()):
            meta = KAL_RUN_META.get(run_name, (None, None, None, None))
            for enst in enst_ids:
                krow = run_data.get(enst, {})
                kal_detail.append({
                    'gene_symbol':   gsym,
                    'enst_id':       enst,
                    'run_name':      run_name,
                    'read_len_bp':   meta[0],
                    'mode':          meta[1],
                    'read_source':   meta[2],
                    'reference':     meta[3],
                    'tpm':           krow.get('tpm', ''),
                    # est_counts: raw EM counts, NOT comparable across runs
                    # (different library sizes); included for reference only
                    'est_counts_DO_NOT_COMPARE_ACROSS_RUNS': krow.get('est_counts', ''),
                    'ambiguity_types': ';'.join(cats),
                })

    if missing_genes:
        print(f'  WARNING: {len(missing_genes)} genes not found in MANE FASTA: '
              f'{", ".join(missing_genes)}')

    return classifications, star_detail, kal_detail
