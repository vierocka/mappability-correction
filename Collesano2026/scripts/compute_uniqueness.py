#!/usr/bin/env python3
"""
Worker for the genomic exon-flank mappability method. Called by
run_mappability_master.sh — never run standalone except for debugging.

Three subcommands, covering the two Step-1 variants (base MANE transcripts
vs. a cell line's personalized HLA-A/B/C) plus the shared Step 3 scoring:

  simulate-base    Step 1 for ~19,373 primary-chromosome MANE transcripts.
                    Genes whose GFF locus is NOT on a primary chromosome
                    (e.g. HLA-DRB4, KIR family, C4A/C4B -- confirmed via the
                    NCBI assembly report to sit on alternate-haplotype MHC
                    scaffolds, i.e. real structural variation, not annotation
                    noise) are EXPLICITLY logged to <outdir>/excluded_non_primary_genes.tsv,
                    never silently dropped. This matches nf-core/rnaseq's own
                    documented recommendation (primary assembly only, to avoid
                    "heavy mapping quality penalties" from alt-locus
                    duplication) -- those genes are not reliably quantifiable
                    from short-read RNA-seq against a single linear reference
                    at all, by the same reasoning the flagship pipeline uses.

  simulate-allele  Step 1 for one cell line's personalized HLA-A/B/C CDS,
                    spliced onto the canonical transcript's genomic exon
                    coordinates. Handles BOTH heterozygous (allele-suffixed
                    ID) and homozygous (original ID, sequence substituted in
                    place -- e.g. JY) entries uniformly; matching on
                    'gene_symbol:HLA-[ABC]' in the header, same convention
                    build_personal_MANE.py itself uses.

  score            Step 3, shared by both modes: count simulated reads per
                    transcript from the R1 FASTQ, scan the aligned BAM for
                    MAPQ=255 (unique) recoveries, write uniqueness_factor.

Genome: ALWAYS the Ensembl-style primary assembly FASTA (bare chromosome
names "1".."22","X","Y","MT", 194 sequences: 25 primary + 169 unplaced/
unlocalized scaffolds, NO alternate haplotypes or patches) -- confirmed
identical between the local copy and the one on HPC
(Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa). The MANE GFF uses
"chr1".."chr22" naming; the chr_map here strips the "chr" prefix, which is
sufficient and unambiguous for this specific file (no accession-mapping
tables needed, unlike the NCBI RefSeq genome).
"""

import argparse
import re
from pathlib import Path
from collections import defaultdict

import pysam
import pandas as pd

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(RC)[::-1]


def parse_fasta(path):
    records = []
    header = seq_parts = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def build_chr_map(genome_fa):
    """chr1 -> '1', chrX -> 'X', chrM/chrMT -> 'MT'. Verified against the
    .fai, not assumed -- refuses to map a chromosome the FASTA doesn't
    actually contain."""
    fai_path = Path(str(genome_fa) + '.fai')
    if not fai_path.exists():
        raise SystemExit(f"Missing .fai: {fai_path} (run samtools faidx first)")
    names = set()
    with open(fai_path) as fh:
        for line in fh:
            names.add(line.split('\t')[0])

    chr_map = {}
    for n in list(range(1, 23)) + ['X', 'Y']:
        if str(n) in names:
            chr_map[f"chr{n}"] = str(n)
    if 'MT' in names:
        chr_map['chrM'] = 'MT'
        chr_map['chrMT'] = 'MT'

    if len(chr_map) < 25:
        raise SystemExit(
            f"Only resolved {len(chr_map)}/25 primary chromosomes from {fai_path} "
            f"-- genome FASTA does not look like the expected primary assembly."
        )
    return chr_map


def gene_name_from_attrs(attrs):
    for part in attrs.split(';'):
        if part.startswith('gene_name='):
            return part.split('=', 1)[1].strip()
    return ''


# ─────────────────────────────────────────────────────────────────────────
# simulate-base
# ─────────────────────────────────────────────────────────────────────────
def cmd_simulate_base(args):
    outdir = Path(args.outdir)
    (outdir / 'simreads').mkdir(parents=True, exist_ok=True)

    chr_map = build_chr_map(args.genome_fa)
    print(f"Primary chromosomes resolved: {len(chr_map)}")

    fasta = pysam.FastaFile(args.genome_fa)

    print("Parsing GFF exon features...")
    tx_exons = defaultdict(list)       # tid -> [(chrom_gff, start0, end, strand)]
    tx_gene_name = {}
    tx_chrom_gff = {}
    with open(args.gff) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            tid = None
            for part in f[8].split(';'):
                if part.startswith('transcript_id='):
                    tid = part.split('=', 1)[1].strip()
                    break
            if tid is None:
                continue
            if f[2] == 'transcript':
                tx_gene_name[tid] = gene_name_from_attrs(f[8])
                tx_chrom_gff[tid] = f[0]
            elif f[2] == 'exon':
                tx_exons[tid].append((f[0], int(f[3]) - 1, int(f[4]), f[6]))

    all_tids = set(tx_exons.keys())
    excluded = []   # (tid, gene_name, chrom_gff) -- not on a primary chromosome
    included_tids = []
    for tid in sorted(all_tids):
        chrom_gff = tx_chrom_gff.get(tid, tx_exons[tid][0][0])
        if chrom_gff in chr_map:
            included_tids.append(tid)
        else:
            excluded.append((tid, tx_gene_name.get(tid, ''), chrom_gff))

    print(f"Transcripts total:            {len(all_tids):,}")
    print(f"On a primary chromosome:      {len(included_tids):,}")
    print(f"Excluded (non-primary locus): {len(excluded):,}")

    excl_path = outdir / 'excluded_non_primary_genes.tsv'
    with open(excl_path, 'w') as fh:
        fh.write("transcript_id\tgene_name\tgff_chromosome\n")
        for tid, gname, chrom in excluded:
            fh.write(f"{tid}\t{gname}\t{chrom}\n")
    print(f"Excluded-gene list written: {excl_path}")

    read_len = args.read_len
    flank = args.flank
    qual = 'I' * read_len
    r1_path = outdir / 'simreads' / 'sim_R1.fastq'
    r2_path = outdir / 'simreads' / 'sim_R2.fastq'
    n_reads = n_skip = 0

    print(f"Simulating reads (READ_LEN={read_len}, FLANK={flank})...")
    with open(r1_path, 'w') as r1, open(r2_path, 'w') as r2:
        for tid in included_tids:
            for ex_idx, (chrom_gff, ex_start, ex_end, strand) in enumerate(tx_exons[tid]):
                chrom_fa = chr_map[chrom_gff]
                win_start = max(0, ex_start - flank)
                win_end = ex_end + flank
                try:
                    seq = fasta.fetch(chrom_fa, win_start, win_end).upper()
                except KeyError:
                    n_skip += 1
                    continue
                L = len(seq)
                if L < read_len:
                    continue
                for pos in range(L - read_len + 1):
                    fwd = seq[pos:pos + read_len]
                    if fwd.count('N') > read_len // 5:
                        continue
                    name = f"{tid}|{ex_idx}|{pos}"
                    r1.write(f"@{name}\n{fwd}\n+\n{qual}\n")
                    r2.write(f"@{name}\n{revcomp(fwd)}\n+\n{qual}\n")
                    n_reads += 1
    fasta.close()

    print(f"Windows skipped (fetch error): {n_skip}")
    print(f"Total reads simulated: {n_reads:,}")
    print(f"Wrote: {r1_path}\nWrote: {r2_path}")


# ─────────────────────────────────────────────────────────────────────────
# simulate-allele
# ─────────────────────────────────────────────────────────────────────────
KNOWN_BASE_IDS = {
    "ENST00000376809.10": "HLA-A",
    "ENST00000412585.7": "HLA-B",
    "ENST00000376228.10": "HLA-C",
}
HLA_RE = re.compile(r'gene_symbol:(HLA-[ABC])\b')


def load_cds_exons(gff_path, transcript_ids):
    """{tid: (chrom_gff, strand, [(start0,end),...] in transcript 5'->3' order)}"""
    exons = defaultdict(list)
    strand_of, chrom_of = {}, {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'CDS':
                continue
            tid = None
            for part in f[8].split(';'):
                if part.startswith('transcript_id='):
                    tid = part.split('=', 1)[1].strip()
                    break
            if tid not in transcript_ids:
                continue
            strand_of[tid] = f[6]
            chrom_of[tid] = f[0]
            exons[tid].append((int(f[3]) - 1, int(f[4])))

    result = {}
    for tid, ex_list in exons.items():
        strand = strand_of[tid]
        ex_sorted = sorted(ex_list, key=lambda x: x[0])
        tx_order = ex_sorted if strand == '+' else list(reversed(ex_sorted))
        result[tid] = (chrom_of[tid], strand, tx_order)
    return result


def cmd_simulate_allele(args):
    outdir = Path(args.outdir)
    (outdir / 'simreads').mkdir(parents=True, exist_ok=True)

    alleles = {}   # id_in_this_sample -> (base_id, seq)
    for header, seq in parse_fasta(args.personal_mane):
        m = HLA_RE.search(header)
        if not m:
            continue
        allele_id = header.split()[0]
        base_id = allele_id if allele_id in KNOWN_BASE_IDS else allele_id.rsplit('_', 1)[0]
        if base_id not in KNOWN_BASE_IDS:
            raise SystemExit(f"Unrecognized base transcript for {allele_id!r}: {base_id!r}")
        alleles[allele_id] = (base_id, seq)

    if not alleles:
        raise SystemExit(f"No HLA-A/B/C records found in {args.personal_mane}")
    print(f"Found {len(alleles)} HLA entry(ies) in {args.personal_mane}:")
    for aid, (base, seq) in alleles.items():
        print(f"  {aid} (base {base}, {len(seq)} bp)")

    cds_exons = load_cds_exons(args.gff, set(KNOWN_BASE_IDS))
    for aid, (base, seq) in alleles.items():
        exon_len_sum = sum(e - s for s, e in cds_exons[base][2])
        if exon_len_sum != len(seq):
            raise SystemExit(
                f"Length mismatch for {aid}: allele CDS {len(seq)} bp vs "
                f"{base} CDS-exon sum {exon_len_sum} bp -- splicing plan invalid."
            )
        print(f"  OK {aid}: {len(seq)} bp matches {base} CDS-exon sum")

    chr_map = build_chr_map(args.genome_fa)
    fasta = pysam.FastaFile(args.genome_fa)

    read_len, flank = args.read_len, args.flank
    qual = 'I' * read_len
    r1_path = outdir / 'simreads' / 'sim_R1.fastq'
    r2_path = outdir / 'simreads' / 'sim_R2.fastq'
    n_reads = n_skip = 0

    with open(r1_path, 'w') as r1, open(r2_path, 'w') as r2:
        for allele_id, (base, seq) in alleles.items():
            chrom_gff, strand, tx_exons = cds_exons[base]
            chrom_fa = chr_map[chrom_gff]

            cursor = 0
            exon_chunks = []
            for (g_start, g_end) in tx_exons:
                elen = g_end - g_start
                exon_chunks.append(seq[cursor:cursor + elen])
                cursor += elen
            assert cursor == len(seq)

            for ex_idx, ((g_start, g_end), chunk) in enumerate(zip(tx_exons, exon_chunks)):
                win_start = max(0, g_start - flank)
                win_end = g_end + flank
                try:
                    geno_seq = fasta.fetch(chrom_fa, win_start, win_end).upper()
                except KeyError:
                    n_skip += 1
                    continue

                left_flank_len = g_start - win_start
                exon_len = g_end - g_start
                left_flank = geno_seq[:left_flank_len]
                right_flank = geno_seq[left_flank_len + exon_len:]
                exon_seq_plus = chunk if strand == '+' else revcomp(chunk)
                spliced = left_flank + exon_seq_plus + right_flank

                L = len(spliced)
                if L < read_len:
                    continue
                for pos in range(L - read_len + 1):
                    fwd = spliced[pos:pos + read_len]
                    if fwd.count('N') > read_len // 5:
                        continue
                    name = f"{allele_id}|{ex_idx}|{pos}"
                    r1.write(f"@{name}\n{fwd}\n+\n{qual}\n")
                    r2.write(f"@{name}\n{revcomp(fwd)}\n+\n{qual}\n")
                    n_reads += 1
    fasta.close()

    print(f"Windows skipped: {n_skip}")
    print(f"Total reads simulated: {n_reads:,}")
    print(f"Wrote: {r1_path}\nWrote: {r2_path}")


# ─────────────────────────────────────────────────────────────────────────
# score  (shared Step 3 -- identical for base and allele modes)
# ─────────────────────────────────────────────────────────────────────────
def cmd_score(args):
    outdir = Path(args.outdir)
    bam_path = Path(args.bam)
    r1_path = outdir / 'simreads' / 'sim_R1.fastq'

    print("Counting simulated reads per transcript...")
    sim_counts = defaultdict(int)
    with open(r1_path) as fh:
        for line in fh:
            if line.startswith('@'):
                sim_counts[line[1:].split('|')[0]] += 1
    print(f"Transcripts: {len(sim_counts):,}  Reads: {sum(sim_counts.values()):,}")

    print("Scanning BAM for MAPQ=255 recoveries...")
    unique_back, total_back = defaultdict(int), defaultdict(int)
    bam = pysam.AlignmentFile(str(bam_path), 'rb')
    for read in bam.fetch(until_eof=True):
        if read.is_supplementary or read.is_secondary or read.is_read2:
            continue
        if read.is_unmapped:
            continue
        tid = read.query_name.split('|')[0]
        total_back[tid] += 1
        if read.mapping_quality == 255:
            unique_back[tid] += 1
    bam.close()

    records = []
    for tid, n_sim in sim_counts.items():
        n_u = unique_back.get(tid, 0)
        n_t = total_back.get(tid, 0)
        records.append({
            'transcript_id': tid,
            'n_positions': n_sim,
            'n_unique_back': n_u,
            'n_multi_back': n_t - n_u,
            'n_unmapped': n_sim - n_t,
            'uniqueness_factor': round(n_u / n_sim, 6) if n_sim > 0 else 0.0,
        })

    df = pd.DataFrame(records).sort_values('uniqueness_factor')
    out = outdir / args.out_name
    df.to_csv(out, sep='\t', index=False)

    print(f"\n── Summary ──────────────────────────────────────")
    print(f"Transcripts scored: {len(df):,}")
    print(f"f=1.00: {(df.uniqueness_factor == 1.0).sum():,}   "
          f"f>=0.90: {(df.uniqueness_factor >= 0.90).sum():,}   "
          f"f<0.10: {(df.uniqueness_factor < 0.10).sum():,}")
    print(f"Saved: {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    p1 = sub.add_parser('simulate-base')
    p1.add_argument('--gff', required=True)
    p1.add_argument('--genome-fa', required=True)
    p1.add_argument('--read-len', type=int, required=True)
    p1.add_argument('--flank', type=int, default=50)
    p1.add_argument('--outdir', required=True)
    p1.set_defaults(func=cmd_simulate_base)

    p2 = sub.add_parser('simulate-allele')
    p2.add_argument('--gff', required=True)
    p2.add_argument('--genome-fa', required=True)
    p2.add_argument('--personal-mane', required=True)
    p2.add_argument('--read-len', type=int, required=True)
    p2.add_argument('--flank', type=int, default=50)
    p2.add_argument('--outdir', required=True)
    p2.set_defaults(func=cmd_simulate_allele)

    p3 = sub.add_parser('score')
    p3.add_argument('--outdir', required=True)
    p3.add_argument('--bam', required=True)
    p3.add_argument('--out-name', default='uniqueness_factors.tsv')
    p3.set_defaults(func=cmd_score)

    args = ap.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
