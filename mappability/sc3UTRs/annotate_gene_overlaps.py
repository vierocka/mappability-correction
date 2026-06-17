#!/usr/bin/env python3
"""
annotate_gene_overlaps.py  —  Genome-wide gene overlap and duplication annotation

For each gene in GENCODE v49, computes:

  1. Genomic overlap group   — connected components of genes sharing genomic coordinates
                                (interval sweep algorithm, per chromosome)

  2. Sequence cluster        — from 3end_clusters_250bp.tsv (transcripts sharing 3' end
                                sequence at >= 90% identity, 250 bp window)

  3. Local duplication group — cluster members on the same chromosome within LOCAL_MB

  4. Phenomenon label        — primary reason for non-uniqueness:

       genomic_overlap_antisense  genes overlapping on opposite strands
       genomic_overlap_sense      genes overlapping on the same strand
       nested                     gene interval fully contained inside another gene
       readthrough                GENCODE readthrough_gene tag
       retrogene                  GENCODE retrogene tag
       local_tandem_dup           same 3' seq cluster, same chrom, within LOCAL_MB
       pseudogene_shadow          PC gene in same seq cluster as a pseudogene, within LOCAL_MB
       distant_paralog            same 3' seq cluster but > LOCAL_MB or different chrom
       unique                     no overlap, no cluster sharing

Inputs:
    GTF_GZ   GENCODE v49 comprehensive annotation
    CLUSTERS  results/gencode_v49/3end_clusters_250bp.tsv
    SUMMARY   results/gencode_v49/3end_cluster_summary_250bp.tsv

Output:
    results/gene_overlap_annotation.tsv

Usage:
    python annotate_gene_overlaps.py
"""

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

# ── paths ──────────────────────────────────────────────────────────────────────
HERE     = Path(__file__).resolve().parent
GTF_GZ   = Path('/home/veve/Dropbox/Self-Nonself/Reference/GENCODE/gencode.v49.annotation.gtf.gz')
CLUSTERS = HERE / 'results/gencode_v49/3end_clusters_250bp.tsv'
SUMMARY  = HERE / 'results/gencode_v49/3end_cluster_summary_250bp.tsv'
OUT      = HERE / 'results/gene_overlap_annotation.tsv'

LOCAL_MB = 5.0     # distance threshold for "local" duplication (Mb)
WINDOW   = 250     # 3' window used for sequence clusters

PSEUDOGENE_TYPES = {
    'processed_pseudogene', 'unprocessed_pseudogene',
    'transcribed_unprocessed_pseudogene', 'transcribed_processed_pseudogene',
    'rRNA_pseudogene', 'transcribed_unitary_pseudogene', 'unitary_pseudogene',
    'IG_V_pseudogene', 'IG_C_pseudogene', 'IG_D_pseudogene', 'IG_J_pseudogene',
    'TR_V_pseudogene', 'TR_J_pseudogene', 'polymorphic_pseudogene',
    'pseudogene',
}


# ── union-find ─────────────────────────────────────────────────────────────────
class UF:
    def __init__(self, n):
        self.p = list(range(n))
        self.r = [0] * n

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return
        if self.r[rx] < self.r[ry]:
            rx, ry = ry, rx
        self.p[ry] = rx
        if self.r[rx] == self.r[ry]:
            self.r[rx] += 1


# ── GTF parsing ────────────────────────────────────────────────────────────────
_re_attr = re.compile(r'(\w+)\s+"([^"]+)"')

def parse_gtf_genes(path):
    genes = []
    open_fn = gzip.open if str(path).endswith('.gz') else open
    with open_fn(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            attr_str = fields[8]
            attrs    = dict(_re_attr.findall(attr_str))
            tags     = re.findall(r'tag "([^"]+)"', attr_str)
            genes.append({
                'gene_id':   attrs.get('gene_id', '').split('.')[0],
                'gene_name': attrs.get('gene_name', ''),
                'gene_type': attrs.get('gene_type', ''),
                'chrom':     fields[0],
                'start':     int(fields[3]),
                'end':       int(fields[4]),
                'strand':    fields[6],
                'tags':      ';'.join(sorted(set(tags))),
            })
    return genes


# ── interval overlap → connected components (per chromosome) ──────────────────
def compute_overlap_groups(genes):
    """
    Returns a dict gene_id → (group_id, group_members_list).
    group_id is the gene_id of the root of the connected component.
    """
    chrom_genes = defaultdict(list)
    for i, g in enumerate(genes):
        chrom_genes[g['chrom']].append(i)

    uf = UF(len(genes))

    for chrom, idxs in chrom_genes.items():
        # sort by start
        sorted_idx = sorted(idxs, key=lambda i: genes[i]['start'])
        # sweep: maintain max_end of the active window
        active = []   # list of (end, idx)
        for i in sorted_idx:
            s = genes[i]['start']
            e = genes[i]['end']
            # prune genes that ended before current start
            active = [(ae, ai) for ae, ai in active if ae >= s]
            for ae, ai in active:
                uf.union(i, ai)
            active.append((e, i))

    # build group membership
    group_members = defaultdict(list)
    for i, g in enumerate(genes):
        group_members[uf.find(i)].append(i)

    return uf, group_members


def overlap_group_stats(gene_indices, genes):
    """Return dict with group statistics."""
    biotypes   = [genes[i]['gene_type'] for i in gene_indices]
    strands    = [genes[i]['strand']    for i in gene_indices]
    names      = [genes[i]['gene_name'] for i in gene_indices]
    biotype_set = set(biotypes)
    strand_set  = set(strands)
    n_pc    = sum(1 for b in biotypes if b == 'protein_coding')
    n_ps    = sum(1 for b in biotypes if b in PSEUDOGENE_TYPES)

    if len(strand_set) == 1:
        strand_cfg = 'same_strand'
    else:
        strand_cfg = 'antisense'

    # check nesting: is any gene fully inside another?
    intervals = [(genes[i]['start'], genes[i]['end'], i) for i in gene_indices]
    intervals.sort()
    nesting_set = set()
    for j, (s1, e1, i1) in enumerate(intervals):
        for s2, e2, i2 in intervals[j+1:]:
            if s2 > e1:
                break
            if s2 >= s1 and e2 <= e1:
                nesting_set.add(i2)
            elif s1 >= s2 and e1 <= e2:
                nesting_set.add(i1)

    return {
        'group_size':            len(gene_indices),
        'group_biotypes':        ';'.join(sorted(biotype_set)),
        'group_strand_config':   strand_cfg,
        'group_n_pc':            n_pc,
        'group_n_pseudogene':    n_ps,
        'group_genes':           ';'.join(sorted(names)),
        'nesting_set':           nesting_set,
    }


# ── 3' sequence clusters ───────────────────────────────────────────────────────
def load_seq_clusters(clusters_path, summary_path, genes_by_id):
    # gene_id → cluster_id (take first cluster seen if multiple isoforms)
    gene_to_cluster = {}
    with open(clusters_path) as fh:
        rdr = csv.DictReader(fh, delimiter='\t')
        for row in rdr:
            gid = row['gene_id'].split('.')[0]
            cid = row['cluster_id']
            if gid not in gene_to_cluster:
                gene_to_cluster[gid] = cid

    # cluster_id → summary stats
    cluster_info = {}
    with open(summary_path) as fh:
        rdr = csv.DictReader(fh, delimiter='\t')
        for row in rdr:
            cluster_info[row['cluster_id']] = {
                'seq_cluster_size':       int(row['cluster_size']),
                'seq_cluster_n_genes':    int(row['n_genes']),
                'seq_cluster_n_pc':       int(row['n_protein_coding']),
                'seq_cluster_n_pseudo':   int(row['n_pseudogene']),
                'seq_cluster_mixed':      int(row['mixed_pc_pseudo']),
                'seq_cluster_gene_names': row['gene_names'],
            }

    return gene_to_cluster, cluster_info


def compute_cluster_genomic_spread(gene_to_cluster, genes):
    """
    For each cluster, gather genomic coords of member genes.
    Returns cluster_id → {chrom: [positions...], multi_chrom: bool, max_span_mb: float}
    """
    cluster_positions = defaultdict(lambda: defaultdict(list))
    for g in genes:
        gid = g['gene_id']
        cid = gene_to_cluster.get(gid)
        if cid:
            cluster_positions[cid][g['chrom']].append(g['start'])
            cluster_positions[cid][g['chrom']].append(g['end'])

    cluster_spread = {}
    for cid, chrom_pos in cluster_positions.items():
        multi_chrom = len(chrom_pos) > 1
        max_span = 0.0
        for chrom, positions in chrom_pos.items():
            span = (max(positions) - min(positions)) / 1e6
            max_span = max(max_span, span)
        cluster_spread[cid] = {
            'multi_chrom': multi_chrom,
            'max_span_mb': round(max_span, 3),
            'n_chroms':    len(chrom_pos),
        }
    return cluster_spread


# ── phenomenon classification ──────────────────────────────────────────────────
def classify_phenomenon(gene, group_stats, cluster_id, cluster_info, cluster_spread, gene_to_cluster, genes_by_id):
    tags     = set(gene['tags'].split(';')) if gene['tags'] else set()
    gtype    = gene['gene_type']
    gsize    = group_stats['group_size']
    strand_c = group_stats['group_strand_config']
    is_nested = (genes_by_id.get(gene['gene_id'], {}).get('idx') in group_stats['nesting_set'])

    # Priority 1: readthrough / retrogene (GENCODE-annotated structural events)
    if 'readthrough_gene' in tags:
        return 'readthrough'
    if 'retrogene' in tags:
        return 'retrogene'

    # Priority 2: genomic overlap
    if gsize > 1:
        if is_nested:
            return 'nested'
        if strand_c == 'antisense':
            return 'genomic_overlap_antisense'
        return 'genomic_overlap_sense'

    # Priority 3: sequence cluster-based (no genomic overlap, but shared sequence)
    if cluster_id and cluster_id in cluster_info:
        ci    = cluster_info[cluster_id]
        cs    = cluster_spread.get(cluster_id, {})
        n_g   = ci['seq_cluster_n_genes']
        n_ps  = ci['seq_cluster_n_pseudo']
        mixed = ci['seq_cluster_mixed']
        multi = cs.get('multi_chrom', True)
        span  = cs.get('max_span_mb', 999.0)

        if n_g > 1:
            if mixed and n_ps > 0 and not multi and span <= LOCAL_MB:
                return 'pseudogene_shadow'
            if not multi and span <= LOCAL_MB:
                return 'local_tandem_dup'
            return 'distant_paralog'

    return 'unique'


# ── main ───────────────────────────────────────────────────────────────────────
def main():
    print('Parsing GENCODE v49 GTF gene features …')
    genes = parse_gtf_genes(GTF_GZ)
    print(f'  {len(genes):,} gene records loaded')

    # index by gene_id (versioned stripped)
    genes_by_id = {}
    for i, g in enumerate(genes):
        g['idx'] = i
        genes_by_id[g['gene_id']] = g

    print('Computing genomic overlap groups …')
    uf, group_members = compute_overlap_groups(genes)

    # precompute per-gene group stats
    print('  Computing group statistics …')
    root_to_stats = {}
    for root, idxs in group_members.items():
        root_to_stats[root] = overlap_group_stats(idxs, genes)

    print('Loading 3\' sequence clusters …')
    gene_to_cluster, cluster_info = load_seq_clusters(CLUSTERS, SUMMARY, genes_by_id)
    print(f'  {len(gene_to_cluster):,} gene→cluster mappings, {len(cluster_info):,} clusters')

    print('Computing cluster genomic spread …')
    cluster_spread = compute_cluster_genomic_spread(gene_to_cluster, genes)

    # ── build output rows ─────────────────────────────────────────────────────
    print('Classifying phenomena and writing output …')
    FIELDS = [
        'gene_id', 'gene_name', 'gene_type',
        'chrom', 'start', 'end', 'strand',
        'gencode_tags',
        # genomic overlap
        'overlap_group_id', 'overlap_group_size',
        'overlap_group_biotypes', 'overlap_strand_config',
        'overlap_n_pc', 'overlap_n_pseudogene',
        'overlap_group_genes',
        # sequence cluster
        'seq_cluster_id',
        'seq_cluster_size', 'seq_cluster_n_genes',
        'seq_cluster_n_pc', 'seq_cluster_n_pseudo', 'seq_cluster_mixed',
        'seq_cluster_multi_chrom', 'seq_cluster_max_span_mb',
        # phenomenon
        'phenomenon',
    ]

    n_written = 0
    with open(OUT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter='\t',
                           lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        for i, g in enumerate(genes):
            root  = uf.find(i)
            gs    = root_to_stats[root]
            cid   = gene_to_cluster.get(g['gene_id'])
            ci    = cluster_info.get(cid, {}) if cid else {}
            cs    = cluster_spread.get(cid, {}) if cid else {}
            ph    = classify_phenomenon(g, gs, cid, cluster_info, cluster_spread, gene_to_cluster, genes_by_id)

            row = {
                'gene_id':               g['gene_id'],
                'gene_name':             g['gene_name'],
                'gene_type':             g['gene_type'],
                'chrom':                 g['chrom'],
                'start':                 g['start'],
                'end':                   g['end'],
                'strand':                g['strand'],
                'gencode_tags':          g['tags'],
                'overlap_group_id':      f'OG_{root}',
                'overlap_group_size':    gs['group_size'],
                'overlap_group_biotypes': gs['group_biotypes'],
                'overlap_strand_config': gs['group_strand_config'],
                'overlap_n_pc':          gs['group_n_pc'],
                'overlap_n_pseudogene':  gs['group_n_pseudogene'],
                'overlap_group_genes':   gs['group_genes'],
                'seq_cluster_id':        cid or '',
                'seq_cluster_size':      ci.get('seq_cluster_size', ''),
                'seq_cluster_n_genes':   ci.get('seq_cluster_n_genes', ''),
                'seq_cluster_n_pc':      ci.get('seq_cluster_n_pc', ''),
                'seq_cluster_n_pseudo':  ci.get('seq_cluster_n_pseudo', ''),
                'seq_cluster_mixed':     ci.get('seq_cluster_mixed', ''),
                'seq_cluster_multi_chrom': cs.get('multi_chrom', ''),
                'seq_cluster_max_span_mb': cs.get('max_span_mb', ''),
                'phenomenon':            ph,
            }
            w.writerow(row)
            n_written += 1

    print(f'  {n_written:,} rows written → {OUT.name}')

    # ── summary ───────────────────────────────────────────────────────────────
    from collections import Counter
    print('\nPhenomenon distribution (all 78k genes):')
    ph_counts = Counter()
    with open(OUT) as fh:
        rdr = csv.DictReader(fh, delimiter='\t')
        for row in rdr:
            ph_counts[row['phenomenon']] += 1
    for ph, n in ph_counts.most_common():
        print(f'  {ph:<35s} {n:>7,}')
    print(f'\nDone → {OUT}')


if __name__ == '__main__':
    main()
