#!/usr/bin/env python3
"""
run_hla_mafft.py — Align HLA transcripts vs their pseudogene relatives using MAFFT.

For each paralog group:
  1. Extracts MANE Select transcript for functional HLAs (longest protein_coding otherwise).
  2. Extracts longest available transcript for each pseudogene.
  3. Writes per-group multi-FASTA.
  4. Runs MAFFT --auto (FASTA + Clustal output).

Source FASTA : GENCODE v49 comprehensive transcripts
               (includes transcribed_unprocessed_pseudogene biotypes)
MANE Select  : MANE.GRCh38.v1.5.ensembl_rna.fna (for functional HLAs)

Output (all in this directory):
  {group}.input.fa      — unaligned sequences
  {group}.aln.fa        — MAFFT alignment, FASTA format
  {group}.aln.clustal   — MAFFT alignment, Clustal format (human-readable)
  {group}.aln.log       — MAFFT stderr log
  transcript_table.tsv  — which transcript was selected for each gene
"""

import gzip, re, shutil, subprocess, sys
from collections import defaultdict
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE     = Path(__file__).resolve().parent
GENCODE  = Path('~/Dropbox/Self-Nonself/Reference/GENCODE/gencode.v49.transcripts.fa.gz').expanduser()
MANE_RNA = Path('~/Dropbox/Self-Nonself/Reference/MANE/MANE.GRCh38.v1.5.ensembl_rna.fna').expanduser()

# ── Paralog groups ────────────────────────────────────────────────────────────
# Each group is aligned as a multiple alignment.
# Pseudogene members are marked with (P) in the alignment header.
GROUPS = {

    # Class I — A clade (HLA-A is the primary reference; H/J/L are its closest relatives)
    'classI_A_clade': {
        'functional':  ['HLA-A', 'HLA-E', 'HLA-G'],
        'pseudogene':  ['HLA-H', 'HLA-J', 'HLA-L', 'HLA-P', 'HLA-V'],
    },

    # Class I — B/C clade
    'classI_BC_clade': {
        'functional':  ['HLA-B', 'HLA-C', 'HLA-F'],
        'pseudogene':  ['HLA-K', 'HLA-N', 'HLA-T', 'HLA-U', 'HLA-W'],
    },

    # Class II — DRB family (all paralogs aligned together)
    'classII_DRB': {
        'functional':  ['HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5'],
        'pseudogene':  ['HLA-DRB2', 'HLA-DRB6', 'HLA-DRB7', 'HLA-DRB8', 'HLA-DRB9'],
    },

    # Class II — DQ locus
    'classII_DQ': {
        'functional':  ['HLA-DQA1', 'HLA-DQB1', 'HLA-DQA2', 'HLA-DQB2'],
        'pseudogene':  ['HLA-DQB3'],
    },

    # Class II — DP locus
    'classII_DP': {
        'functional':  ['HLA-DPA1', 'HLA-DPB1'],
        'pseudogene':  ['HLA-DPA2', 'HLA-DPB2'],
    },
}

# Biotypes that count as pseudogene in GENCODE
PSEUDO_BIOTYPES = {
    'pseudogene', 'unprocessed_pseudogene', 'processed_pseudogene',
    'transcribed_unprocessed_pseudogene', 'transcribed_processed_pseudogene',
    'polymorphic_pseudogene',
}

# ── Load MANE Select transcript IDs ──────────────────────────────────────────

def load_mane_ids(mane_rna_fna):
    """Return dict: gene_symbol → transcript_id (without version) from MANE RNA FASTA."""
    mane = {}
    if not mane_rna_fna.exists():
        print(f'WARNING: MANE RNA not found at {mane_rna_fna}', file=sys.stderr)
        return mane
    with open(mane_rna_fna) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            m = re.search(r'>(ENST\d+)\.\d+.*gene_symbol:(\S+)', line)
            if m:
                tid, gene = m.group(1), m.group(2)
                mane[gene] = tid
    print(f'Loaded {len(mane)} MANE Select transcript IDs')
    return mane


# ── Parse GENCODE transcript FASTA ───────────────────────────────────────────

def parse_gencode(fa_gz, target_genes):
    """
    Read GENCODE transcript FASTA, return:
      gene → list of (transcript_id, transcript_name, length, biotype, seq)
    Only keeps transcripts for genes in target_genes.
    """
    records = defaultdict(list)
    target  = set(target_genes)
    header  = None
    seq_buf = []

    def _flush(hdr, buf):
        if hdr is None:
            return
        parts = hdr.lstrip('>').split('|')
        if len(parts) < 8:
            return
        tid, gene, length_str, biotype = parts[0], parts[5], parts[6], parts[7]
        if gene not in target:
            return
        try:
            length = int(length_str)
        except ValueError:
            length = len(''.join(buf))
        seq = ''.join(buf)
        tname = parts[4] if len(parts) > 4 else tid
        records[gene].append((tid, tname, length, biotype, seq))

    open_fn = gzip.open if str(fa_gz).endswith('.gz') else open
    with open_fn(fa_gz, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                _flush(header, seq_buf)
                header  = line
                seq_buf = []
            else:
                seq_buf.append(line)
    _flush(header, seq_buf)

    print(f'Extracted transcripts for {len(records)} genes from GENCODE')
    return records


# ── Transcript selection ──────────────────────────────────────────────────────

def select_transcript(gene, transcripts, mane_ids, is_pseudogene):
    """
    Functional genes: MANE Select transcript (matched by ENST ID prefix).
    Pseudogenes:      longest transcript (any biotype from PSEUDO_BIOTYPES).
    Fallback:         longest protein_coding or any transcript.
    Returns (tid, tname, length, biotype, seq) or None.
    """
    if not transcripts:
        return None

    if not is_pseudogene and gene in mane_ids:
        mane_id = mane_ids[gene]
        for t in transcripts:
            if t[0].split('.')[0] == mane_id:
                return t
        # MANE transcript not in this FASTA — fall through to longest pc
        print(f'  NOTE: MANE Select {mane_id} for {gene} not in GENCODE FASTA; using longest')

    if is_pseudogene:
        # Prefer transcribed pseudogene biotypes, then any pseudogene
        pseudo = [t for t in transcripts if t[3] in PSEUDO_BIOTYPES]
        if pseudo:
            return max(pseudo, key=lambda t: t[2])

    # Longest protein-coding
    pc = [t for t in transcripts if t[3] == 'protein_coding']
    if pc:
        return max(pc, key=lambda t: t[2])

    # Absolute fallback: longest of anything
    return max(transcripts, key=lambda t: t[2])


# ── FASTA I/O ────────────────────────────────────────────────────────────────

def write_fasta(records, out_path):
    """Write list of (header, seq) to FASTA file, 80 chars per line."""
    with open(out_path, 'w') as fh:
        for header, seq in records:
            fh.write(f'>{header}\n')
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + '\n')


# ── Run MAFFT ────────────────────────────────────────────────────────────────

def run_mafft(input_fa, out_aln_fa, out_aln_clustal, log_path):
    """
    Run MAFFT --auto in FASTA output mode, then re-run for Clustal output.
    For ≤10 sequences we use --localpair --maxiterate 100 (L-INS-i) — more
    accurate for sequences with one conserved domain (HLA peptide-binding groove)
    flanked by divergent UTRs and partially deleted pseudogene exons.
    """
    mafft = shutil.which('mafft')
    if not mafft:
        print('ERROR: mafft not found in PATH', file=sys.stderr)
        sys.exit(1)

    # Count sequences to choose strategy
    n_seqs = sum(1 for l in open(input_fa) if l.startswith('>'))
    strategy = ['--localpair', '--maxiterate', '100'] if n_seqs <= 15 else ['--auto']

    # FASTA output
    with open(log_path, 'w') as log_fh:
        with open(out_aln_fa, 'w') as fa_fh:
            subprocess.run(
                [mafft] + strategy + ['--thread', '-1', str(input_fa)],
                stdout=fa_fh, stderr=log_fh, check=True
            )

    # Clustal output (for human reading — shows conservation track)
    with open(out_aln_clustal, 'w') as cl_fh:
        with open(log_path, 'a') as log_fh:
            subprocess.run(
                [mafft] + strategy + ['--thread', '-1',
                 '--clustalout', str(input_fa)],
                stdout=cl_fh, stderr=log_fh, check=True
            )

    print(f'    → {out_aln_fa.name}  {out_aln_clustal.name}')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    # Collect all target genes across all groups
    all_genes = set()
    pseudo_set = set()
    for grp in GROUPS.values():
        all_genes.update(grp['functional'])
        all_genes.update(grp['pseudogene'])
        pseudo_set.update(grp['pseudogene'])

    mane_ids    = load_mane_ids(MANE_RNA)
    gene_seqs   = parse_gencode(GENCODE, all_genes)

    # Build transcript table
    table_rows = []
    selected   = {}   # gene → (tid, tname, length, biotype, seq)

    for gene in sorted(all_genes):
        is_pseudo = gene in pseudo_set
        transcripts = gene_seqs.get(gene, [])
        t = select_transcript(gene, transcripts, mane_ids, is_pseudo)
        if t is None:
            print(f'  WARNING: no transcript found for {gene}')
        selected[gene] = t
        table_rows.append({
            'gene':        gene,
            'is_pseudogene': is_pseudo,
            'transcript_id': t[0] if t else 'NOT_FOUND',
            'transcript_name': t[1] if t else '',
            'length':      t[2] if t else 0,
            'biotype':     t[3] if t else '',
            'n_isoforms':  len(transcripts),
            'selection_rule': (
                'MANE_Select' if (not is_pseudo and gene in mane_ids and t
                                  and t[0].split('.')[0] == mane_ids.get(gene))
                else 'longest_pseudogene' if is_pseudo
                else 'longest_protein_coding'
            ),
        })

    # Write transcript table
    import csv
    table_path = BASE / 'transcript_table.tsv'
    with open(table_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=[
            'gene', 'is_pseudogene', 'transcript_id', 'transcript_name',
            'length', 'biotype', 'n_isoforms', 'selection_rule'
        ], delimiter='\t')
        w.writeheader()
        w.writerows(table_rows)
    print(f'\nTranscript selection written → {table_path.name}')

    # Per-group alignment
    for group_name, members in GROUPS.items():
        print(f'\n{"═"*60}')
        print(f'Group: {group_name}')
        all_members = members['functional'] + members['pseudogene']
        fasta_records = []
        missing = []

        for gene in all_members:
            t = selected.get(gene)
            if t is None:
                missing.append(gene)
                continue
            is_pseudo = gene in pseudo_set
            tag = '(P)' if is_pseudo else '(F)'
            header = f'{gene} {tag} | {t[0]} | {t[2]}bp | {t[3]}'
            fasta_records.append((header, t[4]))
            print(f'  {gene:<15} {tag}  {t[2]:>5} bp  {t[3]:<40}  {t[0]}')

        if missing:
            print(f'  SKIPPED (no transcript): {", ".join(missing)}')

        if len(fasta_records) < 2:
            print(f'  Fewer than 2 sequences — skipping alignment')
            continue

        input_fa      = BASE / f'{group_name}.input.fa'
        out_aln_fa    = BASE / f'{group_name}.aln.fa'
        out_aln_clust = BASE / f'{group_name}.aln.clustal'
        log_path      = BASE / f'{group_name}.aln.log'

        write_fasta(fasta_records, input_fa)
        run_mafft(input_fa, out_aln_fa, out_aln_clust, log_path)

    print('\nDone.')


if __name__ == '__main__':
    main()
