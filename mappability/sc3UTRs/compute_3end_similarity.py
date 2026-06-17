#!/usr/bin/env python3
"""
3' end sequence similarity across all transcript isoforms.

Three analyses, each run for windows W = 250, 400, 500, 600 bp:

  1. isoform_sharing  — Within each gene: what fraction of isoform pairs share
     an identical last-W bp? Compared to first-W bp (5' window, background).
     If 3' sharing >> 5' sharing, 3' ends are conserved/identical between
     isoforms more than expected from general isoform overlap. This means
     10x reads lose isoform resolution but stay gene-unique. Shared 3' ends
     between isoforms typically arise from alternative 5' start or alternative
     internal splicing leaving the terminal exon unchanged.

  2. pseudogene_xref  — For each protein-coding transcript, does its last-W bp
     exactly match any pseudogene transcript? A match means 10x reads from that
     coding gene can align to the pseudogene locus, reducing apparent expression.
     Outputs the flagged coding gene, the matched pseudogene(s), and window size.

  3. 3end_clusters    — All transcripts (all biotypes) grouped by identical last-W bp.
     Clusters mixing protein_coding and pseudogene biotypes are flagged. Clusters
     that cross gene-family or gene boundaries expose confusion that survives even
     gene-level counting in 10x. Output includes cluster ID, size, biotype
     composition, and gene names.

Two input formats are supported:

  --format gencode (default)
    --fasta  gencode.v49.transcripts.fa.gz   (pipe-delimited GENCODE header)
    --gtf    gencode.v49.basic.annotation.gtf
    IDs: ENST/ENSG. Same namespace as MANE v1.5 in the main pipeline.
    Annotation: GENCODE v49 / Ensembl 115 / GRCh38.p14.
    Exon boundaries may differ slightly from RefSeq models.

  --format refseq  (most internally consistent with the STAR genome index)
    --fasta  transcripts_refseq.fna          (gffread output — plain >NM_/NR_ header)
    --gtf    GRCh38p14_primary.gtf
    IDs: NM_/NR_ (RefSeq). Same annotation source as GCF_000001405.40.
    Exon boundaries are EXACTLY what STAR used for alignment.
    Requires gffread on the cluster (see run_3prime_HPC.sh):
      gffread GRCh38p14_primary.gtf \\
              -g GCF_000001405.40_GRCh38.p14_genomic.fna \\
              -w transcripts_refseq.fna

Usage:
    # GENCODE mode (ready to run — transcripts.fa.gz already on disk):
    python compute_3end_similarity.py \\
        --fasta  /path/to/gencode.v49.transcripts.fa.gz \\
        --gtf    /path/to/gencode.v49.basic.annotation.gtf \\
        --format gencode --label gencode_v49 --outdir results/

    # RefSeq mode (run on cluster after gffread extraction):
    python compute_3end_similarity.py \\
        --fasta  transcripts_refseq.fna \\
        --gtf    /path/to/GRCh38p14_primary.gtf \\
        --format refseq --label refseq_2025 --outdir results/

Output files (one set per window in results/{label}/):
    3end_isoform_sharing_{W}bp.tsv   — per gene, isoform pair sharing rates
    3end_pseudogene_xref_{W}bp.tsv   — coding transcripts matching a pseudogene
    3end_clusters_{W}bp.tsv          — per-transcript cluster membership
    3end_cluster_summary_{W}bp.tsv   — per-cluster summary (size, biotypes, genes)
"""

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from hashlib import md5
from itertools import combinations
from pathlib import Path

WINDOWS = [250, 400, 500, 600]

PSEUDOGENE_TYPES = {
    'processed_pseudogene',
    'unprocessed_pseudogene',
    'transcribed_processed_pseudogene',
    'transcribed_unprocessed_pseudogene',
    'transcribed_unitary_pseudogene',
    'unitary_pseudogene',
    'polymorphic_pseudogene',
    'pseudogene',
    'rRNA_pseudogene',
    'IG_V_pseudogene',
    'IG_C_pseudogene',
    'IG_J_pseudogene',
    'TR_V_pseudogene',
    'TR_J_pseudogene',
}


# ── Argument parsing ───────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--fasta',  required=True,
                   help='Transcript FASTA (GENCODE .fa.gz or gffread-extracted RefSeq .fna)')
    p.add_argument('--gtf',    required=True,
                   help='Annotation GTF (GENCODE or RefSeq)')
    p.add_argument('--format', choices=['gencode', 'refseq'], default='gencode',
                   help='Input format: gencode (pipe-delimited headers) or refseq (plain NM_/NR_ headers)')
    p.add_argument('--label',  default='',
                   help='Label for output subdirectory (e.g. gencode_v49, refseq_2025)')
    p.add_argument('--outdir', default='results',
                   help='Base output directory (default: results/)')
    p.add_argument('--windows', nargs='+', type=int, default=WINDOWS,
                   help=f'3\' window sizes in bp (default: {WINDOWS})')
    return p.parse_args()


# ── FASTA parser for GENCODE format ───────────────────────────────────────────
# Header: >ENST00000832824.1|ENSG00000290825.2|...|gene_name|tx_len|gene_type|
# Fields (pipe-split): [0]=transcript_id, [1]=gene_id, [4]=tx_name, [5]=gene_name,
#                      [6]=tx_len, [7]=gene_type

def open_fasta(path):
    p = Path(path)
    if p.suffix == '.gz':
        return gzip.open(str(p), 'rt')
    return open(str(p))


def parse_gencode_fasta(path):
    """Yield (transcript_id, gene_id, gene_name, gene_type, sequence).
    GENCODE header: >ENST...|ENSG...|...|gene_name|tx_len|gene_type|"""
    tid = gid = gname = gtype = None
    seq_chunks = []
    with open_fasta(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                if tid and seq_chunks:
                    yield tid, gid, gname, gtype, ''.join(seq_chunks)
                seq_chunks = []
                header = line[1:]
                parts = header.split('|')
                tid   = parts[0]
                gid   = parts[1] if len(parts) > 1 else ''
                gname = parts[5] if len(parts) > 5 else ''
                gtype = parts[7] if len(parts) > 7 else ''
                if gtype.endswith('|'):
                    gtype = gtype[:-1]
            else:
                seq_chunks.append(line)
    if tid and seq_chunks:
        yield tid, gid, gname, gtype, ''.join(seq_chunks)


def parse_refseq_fasta(path, tx_to_gene):
    """Yield (transcript_id, gene_id, gene_name, gene_type, sequence).
    RefSeq gffread output: plain >NM_000001.4 or >NR_046018.2 headers.
    tx_to_gene: {transcript_id: (gene_symbol, gene_id, gene_biotype)} from parse_refseq_gtf()."""
    tid = None
    seq_chunks = []
    with open_fasta(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                if tid and seq_chunks:
                    gsym, gid, gtype = tx_to_gene.get(tid, ('', tid, ''))
                    yield tid, gid, gsym, gtype, ''.join(seq_chunks)
                seq_chunks = []
                # gffread output: ">NM_000001.4 ..." — take first token
                tid = line[1:].split()[0]
            else:
                seq_chunks.append(line)
    if tid and seq_chunks:
        gsym, gid, gtype = tx_to_gene.get(tid, ('', tid, ''))
        yield tid, gid, gsym, gtype, ''.join(seq_chunks)


# ── GTF parsers ───────────────────────────────────────────────────────────────

_attr_re = re.compile(r'(\w+)\s+"([^"]+)"')


def parse_gtf_gene_types(gtf_path):
    """Return {gene_id_base: gene_type} from gene-level GTF records (GENCODE or RefSeq)."""
    gene_types = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t')
            if len(f) < 9 or f[2] != 'gene':
                continue
            attrs = dict(_attr_re.findall(f[8]))
            gid = attrs.get('gene_id', '').split('.')[0]
            gt  = attrs.get('gene_type', attrs.get('gene_biotype', ''))
            if not gt:
                # RefSeq uses pseudo="true" instead of a biotype field on gene records
                if attrs.get('pseudo') == 'true':
                    gbkey = attrs.get('gbkey', '')
                    gt = 'pseudogene' if gbkey else 'pseudogene'
            if gid:
                gene_types[gid] = gt
    return gene_types


def parse_refseq_gtf(gtf_path):
    """Parse RefSeq GTF; return ({transcript_id: (gene_symbol, gene_id, gene_biotype)},
                                   {gene_id: gene_biotype}).
    RefSeq uses NC_ chromosome names; gene_biotype is on the gene record,
    not on transcript records. Pseudogenes are identified by pseudo="true"."""
    gene_biotypes = {}   # {gene_id: biotype}
    tx_to_gene    = {}   # {transcript_id: (gene_symbol, gene_id, gene_biotype)}

    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t')
            if len(f) < 9:
                continue
            feat  = f[2]
            attrs = dict(_attr_re.findall(f[8]))

            if feat == 'gene':
                gid = attrs.get('gene_id', '')
                gt  = attrs.get('gene_biotype', '')
                if not gt:
                    if attrs.get('pseudo') == 'true':
                        gt = 'pseudogene'
                    else:
                        gbkey = attrs.get('gbkey', '')
                        gt = {'Gene': 'protein_coding', 'misc_RNA': 'ncRNA'}.get(gbkey, gbkey)
                if gid:
                    gene_biotypes[gid] = gt

            elif feat == 'transcript':
                tid = attrs.get('transcript_id', '')
                gid = attrs.get('gene_id', '')
                gsym = attrs.get('gene', gid)
                if tid and gid:
                    tx_to_gene[tid] = (gsym, gid, '')   # biotype filled below

    # Fill in biotype from gene records
    for tid, (gsym, gid, _) in tx_to_gene.items():
        tx_to_gene[tid] = (gsym, gid, gene_biotypes.get(gid, ''))

    return tx_to_gene, gene_biotypes


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    label  = args.label or args.format
    outdir = Path(args.outdir) / label
    outdir.mkdir(parents=True, exist_ok=True)
    windows = sorted(args.windows)

    print(f"\n{'═'*64}")
    print(f"3' end similarity — {label}")
    print(f"  Format : {args.format}")
    print(f"  FASTA  : {args.fasta}")
    print(f"  GTF    : {args.gtf}")
    print(f"  Windows: {windows} bp")
    print(f"  Outdir : {outdir}")
    print(f"{'═'*64}", flush=True)

    # ── Load annotation ────────────────────────────────────────────────────────
    print(f"\nLoading GTF ({args.format} format) ...", flush=True)
    if args.format == 'refseq':
        tx_to_gene, gtf_gene_types = parse_refseq_gtf(args.gtf)
        print(f"  {len(tx_to_gene):,} transcripts, {len(gtf_gene_types):,} genes", flush=True)
    else:
        gtf_gene_types = parse_gtf_gene_types(args.gtf)
        tx_to_gene     = {}
        print(f"  {len(gtf_gene_types):,} gene_ids loaded", flush=True)

    # ── Load all transcript sequences ──────────────────────────────────────────
    print("\nReading transcript FASTA ...", flush=True)
    tx_meta = {}        # {tid: (gid_base, gene_name, gene_type)}
    tx_seq  = {}        # {tid: sequence}
    gene_to_tids = defaultdict(list)

    fasta_iter = (parse_refseq_fasta(args.fasta, tx_to_gene)
                  if args.format == 'refseq'
                  else parse_gencode_fasta(args.fasta))

    n = 0
    for tid, gid, gname, gtype, seq in fasta_iter:
        gid_base = gid.split('.')[0]
        tid_base = tid

        if not gtype and gid_base in gtf_gene_types:
            gtype = gtf_gene_types[gid_base]

        tx_meta[tid_base] = (gid_base, gname, gtype)
        tx_seq[tid_base]  = seq
        gene_to_tids[gid_base].append(tid_base)
        n += 1
        if n % 50_000 == 0:
            print(f"  {n:,} transcripts ...", flush=True)

    print(f"  Total: {n:,} transcripts, {len(gene_to_tids):,} genes", flush=True)

    # ── Per-window analyses ────────────────────────────────────────────────────
    for W in windows:
        print(f"\n{'─'*64}")
        print(f"Window W = {W} bp", flush=True)

        # Build 3' and 5' end sequences per transcript
        ends_3p = {}   # {tid: last-W bp}
        ends_5p = {}   # {tid: first-W bp}
        for tid, seq in tx_seq.items():
            ends_3p[tid] = seq[-W:]
            ends_5p[tid] = seq[:W]

        # ── Analysis 1: isoform sharing ────────────────────────────────────────
        print("  Analysis 1: isoform 3'/5' sharing per gene ...", flush=True)
        iso_rows = []
        for gid_base, tids in gene_to_tids.items():
            if len(tids) < 2:
                continue
            gname = tx_meta[tids[0]][1]
            gtype = tx_meta[tids[0]][2]
            pairs = list(combinations(tids, 2))
            n_pairs = len(pairs)
            n_share_3p = sum(1 for a, b in pairs
                             if ends_3p.get(a) and ends_3p.get(b)
                             and ends_3p[a] == ends_3p[b])
            n_share_5p = sum(1 for a, b in pairs
                             if ends_5p.get(a) and ends_5p.get(b)
                             and ends_5p[a] == ends_5p[b])
            iso_rows.append({
                'gene_id':         gid_base,
                'gene_name':       gname,
                'gene_type':       gtype,
                'n_isoforms':      len(tids),
                'n_pairs':         n_pairs,
                'n_share_3p':      n_share_3p,
                'n_share_5p':      n_share_5p,
                'rate_3p':         round(n_share_3p / n_pairs, 6) if n_pairs else 0.0,
                'rate_5p':         round(n_share_5p / n_pairs, 6) if n_pairs else 0.0,
                'delta_3p_vs_5p':  round((n_share_3p - n_share_5p) / n_pairs, 6) if n_pairs else 0.0,
            })

        iso_rows.sort(key=lambda r: -r['rate_3p'])
        iso_path = outdir / f'3end_isoform_sharing_{W}bp.tsv'
        iso_fields = ['gene_id','gene_name','gene_type','n_isoforms','n_pairs',
                      'n_share_3p','n_share_5p','rate_3p','rate_5p','delta_3p_vs_5p']
        with open(iso_path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=iso_fields, delimiter='\t')
            w.writeheader()
            w.writerows(iso_rows)
        n_genes_full = sum(1 for r in iso_rows if r['rate_3p'] == 1.0)
        n_genes_any  = sum(1 for r in iso_rows if r['n_share_3p'] > 0)
        n_genes_more_3p = sum(1 for r in iso_rows if r['delta_3p_vs_5p'] > 0)
        print(f"    Multi-isoform genes: {len(iso_rows):,}")
        print(f"    All isoforms share 3' end: {n_genes_full:,}")
        print(f"    Any isoform pair shares 3' end: {n_genes_any:,}")
        print(f"    3' sharing > 5' sharing: {n_genes_more_3p:,}")
        print(f"    Saved → {iso_path}", flush=True)

        # ── Analysis 2: pseudogene cross-reference ─────────────────────────────
        print("  Analysis 2: pseudogene 3' end cross-reference ...", flush=True)

        # Build {3'_seq: [(tid, gene_name, gtype), ...]} for pseudogenes
        pseudo_3p = defaultdict(list)
        for tid, seq in tx_seq.items():
            gtype = tx_meta[tid][2]
            if gtype in PSEUDOGENE_TYPES:
                s = ends_3p.get(tid, '')
                if s:
                    pseudo_3p[s].append((tid, tx_meta[tid][1], gtype))

        pseudo_rows = []
        for tid, seq in tx_seq.items():
            gtype = tx_meta[tid][2]
            if gtype != 'protein_coding':
                continue
            s = ends_3p.get(tid, '')
            if not s:
                continue
            matches = pseudo_3p.get(s)
            if not matches:
                continue
            gid_base, gname, _ = tx_meta[tid]
            for m_tid, m_gene, m_type in matches:
                m_gid = tx_meta[m_tid][0]
                pseudo_rows.append({
                    'transcript_id':         tid,
                    'gene_id':               gid_base,
                    'gene_name':             gname,
                    'pseudogene_transcript': m_tid,
                    'pseudogene_gene_id':    m_gid,
                    'pseudogene_gene_name':  m_gene,
                    'pseudogene_type':       m_type,
                    'window_bp':             W,
                    'shared_seq_len':        len(s),
                })

        pseudo_rows.sort(key=lambda r: (r['gene_name'], r['pseudogene_gene_name']))
        pseudo_path = outdir / f'3end_pseudogene_xref_{W}bp.tsv'
        pseudo_fields = ['transcript_id','gene_id','gene_name',
                         'pseudogene_transcript','pseudogene_gene_id','pseudogene_gene_name',
                         'pseudogene_type','window_bp','shared_seq_len']
        with open(pseudo_path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=pseudo_fields, delimiter='\t')
            w.writeheader()
            w.writerows(pseudo_rows)
        n_pc_affected = len({r['gene_id'] for r in pseudo_rows})
        n_pseudo_used = len({r['pseudogene_gene_id'] for r in pseudo_rows})
        print(f"    Protein-coding genes with pseudogene 3' match: {n_pc_affected:,}")
        print(f"    Pseudogenes involved: {n_pseudo_used:,}")
        print(f"    Saved → {pseudo_path}", flush=True)

        # ── Analysis 3: 3' end clusters ────────────────────────────────────────
        print("  Analysis 3: 3' end clusters ...", flush=True)

        # Group all transcripts by identical 3' end sequence
        seq_to_tids = defaultdict(list)
        for tid, s in ends_3p.items():
            if s:
                seq_to_tids[s].append(tid)

        # Keep only clusters with ≥2 members
        clusters = [(s, tids) for s, tids in seq_to_tids.items() if len(tids) > 1]
        # Sort by cluster size descending
        clusters.sort(key=lambda x: -len(x[1]))

        cluster_tx_rows   = []
        cluster_summ_rows = []
        for cid, (seq, tids) in enumerate(clusters, 1):
            seq_hash = md5(seq.encode()).hexdigest()[:10]
            n_tx = len(tids)
            biotype_counts = defaultdict(int)
            gene_ids  = set()
            gene_names = set()
            for tid in tids:
                gid, gname, gtype = tx_meta[tid]
                biotype_counts[gtype] += 1
                gene_ids.add(gid)
                gene_names.add(gname)
            n_pc    = biotype_counts.get('protein_coding', 0)
            n_pseudo = sum(biotype_counts.get(t, 0) for t in PSEUDOGENE_TYPES)
            mixed   = n_pc > 0 and n_pseudo > 0
            for tid in tids:
                gid, gname, gtype = tx_meta[tid]
                cluster_tx_rows.append({
                    'cluster_id':      f'C{cid}_{seq_hash}',
                    'cluster_size':    n_tx,
                    'transcript_id':   tid,
                    'gene_id':         gid,
                    'gene_name':       gname,
                    'gene_type':       gtype,
                    'mixed_pc_pseudo': int(mixed),
                    'window_bp':       W,
                })
            cluster_summ_rows.append({
                'cluster_id':      f'C{cid}_{seq_hash}',
                'cluster_size':    n_tx,
                'n_genes':         len(gene_ids),
                'n_protein_coding': n_pc,
                'n_pseudogene':    n_pseudo,
                'n_other':         n_tx - n_pc - n_pseudo,
                'mixed_pc_pseudo': int(mixed),
                'gene_names':      ';'.join(sorted(gene_names)[:20]),
                'window_bp':       W,
            })

        tx_path   = outdir / f'3end_clusters_{W}bp.tsv'
        summ_path = outdir / f'3end_cluster_summary_{W}bp.tsv'
        tx_fields = ['cluster_id','cluster_size','transcript_id','gene_id',
                     'gene_name','gene_type','mixed_pc_pseudo','window_bp']
        summ_fields = ['cluster_id','cluster_size','n_genes','n_protein_coding',
                       'n_pseudogene','n_other','mixed_pc_pseudo','gene_names','window_bp']
        with open(tx_path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=tx_fields, delimiter='\t')
            w.writeheader()
            w.writerows(cluster_tx_rows)
        with open(summ_path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=summ_fields, delimiter='\t')
            w.writeheader()
            w.writerows(cluster_summ_rows)

        n_clusters     = len(clusters)
        n_mixed        = sum(1 for r in cluster_summ_rows if r['mixed_pc_pseudo'])
        n_cross_gene   = sum(1 for r in cluster_summ_rows if r['n_genes'] > 1)
        n_tx_in_clust  = len(cluster_tx_rows)
        print(f"    Multi-member clusters: {n_clusters:,}")
        print(f"    Clusters with >1 gene: {n_cross_gene:,}")
        print(f"    Clusters mixing PC+pseudogene: {n_mixed:,}")
        print(f"    Transcripts in clusters: {n_tx_in_clust:,}")
        print(f"    Saved → {tx_path}")
        print(f"    Saved → {summ_path}", flush=True)

    print(f"\n{'═'*64}")
    print("All windows done.")
    print(f"Output directory: {outdir.resolve()}")


if __name__ == '__main__':
    main()
