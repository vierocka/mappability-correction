#!/usr/bin/env python3
"""
Generate the 30 missing mappability correction bash scripts.

Each script is built from sections that depend on:
  ref       : genome | dedup_cds | mane
  src       : exon_flank | mane_rna
  mode      : PE | SE
  read_len  : 75 | 100 | 150 | 200
"""

import os, stat
from pathlib import Path

HERE = Path(__file__).resolve().parent

# ── New scripts to generate ────────────────────────────────────────────────────
# (script_id, ref, src, mode, read_len)
NEW_SCRIPTS = [
    # Full genome
    ('v1_d', 'genome',    'mane_rna',   'PE', 100),
    ('v3_b', 'genome',    'exon_flank', 'SE',  75),
    ('v3_c', 'genome',    'mane_rna',   'SE',  75),
    ('v3_d', 'genome',    'mane_rna',   'PE',  75),
    ('v4_b', 'genome',    'exon_flank', 'SE', 150),
    ('v4_c', 'genome',    'mane_rna',   'SE', 150),
    ('v4_d', 'genome',    'mane_rna',   'PE', 150),
    ('v6_a', 'genome',    'mane_rna',   'SE', 200),
    ('v6_b', 'genome',    'mane_rna',   'PE', 200),
    # Dedup CDS
    ('v2_e', 'dedup_cds', 'mane_rna',   'PE', 100),
    ('v2_f', 'dedup_cds', 'exon_flank', 'PE',  75),
    ('v2_g', 'dedup_cds', 'exon_flank', 'SE',  75),
    ('v2_h', 'dedup_cds', 'exon_flank', 'PE', 150),
    ('v2_i', 'dedup_cds', 'exon_flank', 'SE', 150),
    ('v2_j', 'dedup_cds', 'mane_rna',   'SE',  75),
    ('v2_k', 'dedup_cds', 'mane_rna',   'PE',  75),
    ('v2_l', 'dedup_cds', 'mane_rna',   'SE', 150),
    ('v2_m', 'dedup_cds', 'mane_rna',   'PE', 150),
    ('v2_n', 'dedup_cds', 'mane_rna',   'SE', 200),
    ('v2_o', 'dedup_cds', 'mane_rna',   'PE', 200),
    # MANE ref
    ('v5_e', 'mane',      'mane_rna',   'SE', 150),
    ('v5_f', 'mane',      'mane_rna',   'PE',  75),
    ('v5_g', 'mane',      'mane_rna',   'PE', 100),
    ('v5_h', 'mane',      'mane_rna',   'PE', 150),
    ('v5_i', 'mane',      'mane_rna',   'PE', 200),
    ('v5_j', 'mane',      'exon_flank', 'SE',  75),
    ('v5_k', 'mane',      'exon_flank', 'SE', 100),
    ('v5_l', 'mane',      'exon_flank', 'SE', 150),
    ('v5_m', 'mane',      'exon_flank', 'PE',  75),
    ('v5_n', 'mane',      'exon_flank', 'PE', 150),
]

# ── Naming helpers ─────────────────────────────────────────────────────────────
REF_LABELS  = {'genome': 'full genome', 'dedup_cds': 'dedup CDS', 'mane': 'MANE ref'}
SRC_LABELS  = {'exon_flank': 'exon-flank', 'mane_rna': 'MANE RNA'}

def outdir_suffix(ref, src, mode, L):
    """Results subdirectory name (without leading 'mappability_')."""
    L_tag = '' if L == 100 else f'_L{L}bp'
    if ref == 'genome':
        if src == 'exon_flank':
            return f'genomic{"" if mode=="PE" else "_SE"}{L_tag}'
        else:
            return f'genomic_RNA_{mode}{L_tag}'
    elif ref == 'dedup_cds':
        if src == 'exon_flank':
            return f'dedup_cds_{mode}_flank{L_tag}'
        else:
            return f'dedup_cds_RNA_{mode}{L_tag}'
    else:  # mane
        if src == 'exon_flank':
            return f'MANE_exonflank_{mode}{L_tag}'
        else:
            return f'MANE_RNA_{mode}{L_tag}'

def outdir_name(ref, src, mode, L):
    return 'mappability_' + outdir_suffix(ref, src, mode, L)

def tsv_filename(ref, src, mode, L):
    return 'transcript_uniqueness_factors_' + outdir_suffix(ref, src, mode, L) + '.tsv'


# ── Script sections ────────────────────────────────────────────────────────────

def section_header(sid, ref, src, mode, L):
    ref_l = REF_LABELS[ref]; src_l = SRC_LABELS[src]
    return f'''#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — {sid}
#
# Reference : {ref_l}
# Reads     : {src_l}
# Mode      : {mode}
# Read len  : {L} bp
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
REF="$SCRIPT_DIR/../Ref"
MANE_GFF="$REF/MANE.GRCh38.v1.5.ensembl_genomic.gff"
GENOME_FA="$REF/GCF_000001405.40_GRCh38.p14_genomic.fna"
MANE_RNA="$REF/MANE.GRCh38.v1.5.ensembl_rna.fna"
DEDUP_FA="$REF/GCF_000001405.40_GRCh38.p14_cds_from_genomic.dedup.fna"
GENOME_IDX="$REF/STAR_index_genome"
CDS_IDX="$REF/STAR_index_dedup_cds"
MANE_IDX="$REF/STAR_index_MANE"
OUTDIR="$SCRIPT_DIR/results/{outdir_name(ref, src, mode, L)}"
THREADS=64
READ_LEN={L}
FLANK=50

mkdir -p "$OUTDIR/simreads" "$OUTDIR/bam" "$OUTDIR/logs"

echo "════════════════════════════════════════════════════════"
echo "Mappability — {sid} ({ref_l}, {src_l}, {mode}, {L} bp)"
echo "READ_LEN=$READ_LEN  THREADS=$THREADS"
echo "════════════════════════════════════════════════════════"
'''

def section_chr_map():
    """Step 0 — chr→NC map (exon-flank scripts only)."""
    return r'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 0: Build chr→NC chromosome name map
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0: Building chr→NC chromosome name map..."
python3 - << PYEOF
import re
from pathlib import Path

genome_fa = Path("$GENOME_FA")
out_map   = Path("$OUTDIR/chr_name_map.tsv")

ncbi_map = {}
print("  Scanning FASTA headers...", flush=True)
with open(genome_fa) as fh:
    for line in fh:
        if not line.startswith('>'): continue
        name = line[1:].split()[0]
        rest = line.strip()
        m = re.search(r'chromosome (\w+)[,\s]', rest)
        if m and name.startswith('NC_'):
            ncbi_map[f'chr{m.group(1)}'] = name
        if 'mitochondrion' in rest.lower() or 'mitochondrial' in rest.lower():
            ncbi_map['chrM'] = name; ncbi_map['chrMT'] = name

with open(out_map, 'w') as fh:
    for k, v in sorted(ncbi_map.items()):
        fh.write(f"{k}\t{v}\n")
print(f"  Chromosomes mapped: {len(ncbi_map)}")
PYEOF
echo "Step 0 complete."
'''

def section_build_index(ref):
    """Step 0b — build STAR index (skipped if already exists)."""
    if ref == 'genome':
        idx_var, idx_path = 'GENOME_IDX', '"$GENOME_IDX"'
        fa_var = '"$GENOME_FA"'
        sa, bins = 14, 14
        label = 'full genome'
        log = '"$OUTDIR/logs/star_index.log"'
    elif ref == 'dedup_cds':
        idx_var, idx_path = 'CDS_IDX', '"$CDS_IDX"'
        fa_var = '"$DEDUP_FA"'
        sa, bins = 12, 11
        label = 'dedup CDS'
        log = '"$OUTDIR/logs/star_index.log"'
    else:  # mane
        idx_var, idx_path = 'MANE_IDX', '"$MANE_IDX"'
        fa_var = '"$MANE_RNA"'
        sa, bins = 12, 11
        label = 'MANE RNA'
        log = '"$OUTDIR/logs/star_index_MANE.log"'

    return f'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 0b: Build STAR index ({label}) — shared; skipped if already exists
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0b: Checking STAR index ({label})..."
mkdir -p {idx_path}
if [ ! -f {idx_path}/SA ]; then
    STAR \\
        --runMode             genomeGenerate \\
        --genomeDir           {idx_path} \\
        --genomeFastaFiles    {fa_var} \\
        --genomeSAindexNbases {sa} \\
        --genomeChrBinNbits   {bins} \\
        --runThreadN          $THREADS \\
        2>&1 | tee {log}
    echo "Step 0b complete."
else
    echo "Step 0b: index already exists, skipping."
fi
'''

def section_simulate_exon_flank(mode):
    """Step 1 — exon-flank read simulation."""
    if mode == 'PE':
        open_files   = 'with open(r1_path, \'w\') as r1, open(r2_path, \'w\') as r2:'
        r2_decl      = "r2_path = outdir / 'simreads/sim_R2.fastq'"
        write_reads  = (
            "                r1.write(f'@{name}\\n{fwd}\\n+\\n{qual}\\n')\n"
            "                r2.write(f'@{name}\\n{revcomp(fwd)}\\n+\\n{qual}\\n')"
        )
        revcomp_func = "RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')\ndef revcomp(s): return s.translate(RC)[::-1]\n"
    else:
        open_files   = 'with open(r1_path, \'w\') as r1:'
        r2_decl      = ''
        write_reads  = "                r1.write(f'@{name}\\n{fwd}\\n+\\n{qual}\\n')"
        revcomp_func = ''

    return f'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Parse GFF, fetch exon ± FLANK from genome, simulate reads
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Simulating exon-flank reads (READ_LEN=$READ_LEN, FLANK=$FLANK, {mode})..."

python3 - << PYEOF
import subprocess, re
from pathlib import Path
from collections import defaultdict

gff_path  = Path("$MANE_GFF")
genome_fa = Path("$GENOME_FA")
outdir    = Path("$OUTDIR")
read_len  = int("$READ_LEN")
flank     = int("$FLANK")

{revcomp_func}
chr_map = {{}}
map_file = outdir / "chr_name_map.tsv"
with open(map_file) as fh:
    for line in fh:
        parts = line.strip().split('\\t')
        if len(parts) == 2: chr_map[parts[0]] = parts[1]
print(f"  Loaded {{len(chr_map)}} chr->NC mappings", flush=True)

tx_exons = defaultdict(list)
with open(gff_path) as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip('\\n').split('\\t')
        if len(f) < 9 or f[2] != 'exon': continue
        chrom = chr_map.get(f[0], f[0])
        start = int(f[3]) - 1; end = int(f[4]); strand = f[6]
        tid = None
        for part in f[8].split(';'):
            if part.startswith('transcript_id='):
                tid = part.split('=', 1)[1].strip(); break
        if tid: tx_exons[tid].append((chrom, start, end, strand))
print(f"  Transcripts: {{len(tx_exons):,}}", flush=True)

fai = Path(str(genome_fa) + '.fai')
if not fai.exists():
    import subprocess as sp
    sp.run(['samtools', 'faidx', str(genome_fa)], check=True)

r1_path = outdir / 'simreads/sim_R1.fastq'
{r2_decl}
qual = 'I' * read_len
n_reads = n_skip = 0

{open_files}
    for tid, exons in tx_exons.items():
        for ex_idx, (chrom, ex_start, ex_end, strand) in enumerate(exons):
            win_start = max(0, ex_start - flank)
            win_end   = ex_end + flank
            region    = f"{{chrom}}:{{win_start + 1}}-{{win_end}}"
            try:
                res = subprocess.run(['samtools', 'faidx', str(genome_fa), region],
                                     capture_output=True, text=True, check=True)
                seq = ''.join(res.stdout.split('\\n')[1:]).upper()
            except subprocess.CalledProcessError:
                n_skip += 1; continue
            L = len(seq)
            if L < read_len: continue
            for pos in range(L - read_len + 1):
                fwd = seq[pos : pos + read_len]
                if fwd.count('N') > read_len // 5: continue
                name = f"{{tid}}|{{ex_idx}}|{{pos}}"
{write_reads}
                n_reads += 1

print(f"  Exon windows skipped: {{n_skip}}")
print(f"  Total reads simulated: {{n_reads:,}}", flush=True)
PYEOF
echo "Step 1 complete."
'''

def section_simulate_mane_rna(mode):
    """Step 1 — MANE RNA read simulation."""
    if mode == 'PE':
        r2_decl     = "r2_path = outdir / 'simreads/sim_R2.fastq'"
        open_files  = "with open(r1_path, 'w') as r1, open(r2_path, 'w') as r2:"
        write_reads = (
            "            r1.write(f'@{tid}|{pos}\\n{fwd}\\n+\\n{qual}\\n')\n"
            "            r2.write(f'@{tid}|{pos}\\n{revcomp(fwd)}\\n+\\n{qual}\\n')"
        )
        revcomp_func = "RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')\ndef revcomp(s): return s.translate(RC)[::-1]\n"
    else:
        r2_decl      = ''
        open_files   = "with open(r1_path, 'w') as r1:"
        write_reads  = "            r1.write(f'@{tid}|{pos}\\n{fwd}\\n+\\n{qual}\\n')"
        revcomp_func = ''

    return f'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Simulate reads from MANE spliced RNA sequences
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Simulating MANE RNA reads (READ_LEN=$READ_LEN, {mode})..."

python3 - << PYEOF
from pathlib import Path

mane_rna = Path("$MANE_RNA")
outdir   = Path("$OUTDIR")
read_len = int("$READ_LEN")

{revcomp_func}
def parse_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name: yield name, ''.join(seq)
                name = line[1:].split()[0]; seq = []
            else: seq.append(line.upper())
    if name: yield name, ''.join(seq)

r1_path = outdir / 'simreads/sim_R1.fastq'
{r2_decl}
qual = 'I' * read_len
n_reads = n_tx = n_skip = 0

{open_files}
    for tid, seq in parse_fasta(mane_rna):
        if len(seq) < read_len: n_skip += 1; continue
        n_tx += 1
        for pos in range(len(seq) - read_len + 1):
            fwd = seq[pos : pos + read_len]
            if fwd.count('N') > read_len // 5: continue
{write_reads}
            n_reads += 1

print(f"  Transcripts: {{n_tx:,}}  (skipped too short: {{n_skip}})")
print(f"  Total reads simulated: {{n_reads:,}}", flush=True)
PYEOF
echo "Step 1 complete."
'''

def section_star_align(ref, mode, outdir_n):
    """Step 2 — STAR alignment."""
    if ref == 'genome':
        idx = '"$GENOME_IDX"'
        extra = ''
        label = 'full genome'
    elif ref == 'dedup_cds':
        idx = '"$CDS_IDX"'
        extra = ''
        label = 'dedup CDS'
    else:  # mane
        idx = '"$MANE_IDX"'
        extra = '    --alignIntronMax        1 \\\n    --alignEndsType         EndToEnd \\\n'
        label = 'MANE ref'

    if mode == 'PE':
        reads_in = (
            '    --readFilesIn           "$OUTDIR/simreads/sim_R1.fastq" \\\n'
            '                            "$OUTDIR/simreads/sim_R2.fastq" \\'
        )
    else:
        reads_in = '    --readFilesIn           "$OUTDIR/simreads/sim_R1.fastq" \\'

    return f'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: STAR alignment — {mode}, {label}
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 2: STAR alignment ({mode}, {label}, $THREADS threads)..."

STAR \\
    --runThreadN            $THREADS \\
    --genomeDir             {idx} \\
{reads_in}
    --outSAMtype            BAM SortedByCoordinate \\
    --outSAMattributes      NH HI AS NM \\
    --outSAMmultNmax        1 \\
    --outFilterMultimapNmax 40 \\
{extra}    --outBAMsortingThreadN  $THREADS \\
    --outBAMsortingBinsN    20 \\
    --limitBAMsortRAM       160000000000 \\
    --outFileNamePrefix     "$OUTDIR/bam/sim_" \\
    2>&1 | tee "$OUTDIR/logs/star_sim.log"

samtools index -@ $THREADS "$OUTDIR/bam/sim_Aligned.sortedByCoord.out.bam"
echo "Step 2 complete."
'''

def section_compute_uniqueness(mode, tsv_n):
    """Step 3 — compute uniqueness factors."""
    if mode == 'PE':
        bam_filter = 'if read.is_supplementary or read.is_secondary or read.is_read2: continue'
    else:
        bam_filter = 'if read.is_supplementary or read.is_secondary or read.is_unmapped: continue'

    return f'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Compute per-transcript uniqueness factor
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 3: Computing uniqueness factors..."

python3 - << PYEOF
import pysam
from collections import defaultdict
from pathlib import Path

outdir   = Path("$OUTDIR")
bam_path = outdir / "bam/sim_Aligned.sortedByCoord.out.bam"
r1_path  = outdir / "simreads/sim_R1.fastq"

sim_counts = defaultdict(int)
with open(r1_path) as fh:
    for line in fh:
        if line.startswith('@'):
            sim_counts[line[1:].split('|')[0]] += 1
print(f"  Transcripts: {{len(sim_counts):,}}  Reads: {{sum(sim_counts.values()):,}}", flush=True)

unique_back = defaultdict(int)
total_back  = defaultdict(int)
bam = pysam.AlignmentFile(str(bam_path), 'rb')
for read in bam.fetch():
    {bam_filter}
    tid = read.query_name.split('|')[0]
    total_back[tid] += 1
    if read.mapping_quality == 255: unique_back[tid] += 1
bam.close()

import csv
records = []
for tid, n_sim in sim_counts.items():
    n_u = unique_back.get(tid, 0); n_t = total_back.get(tid, 0)
    records.append({{
        'transcript_id':     tid,
        'n_positions':       n_sim,
        'n_unique_back':     n_u,
        'n_multi_back':      n_t - n_u,
        'n_unmapped':        n_sim - n_t,
        'uniqueness_factor': round(n_u / n_sim, 6) if n_sim > 0 else 0.0,
    }})

records.sort(key=lambda r: r['uniqueness_factor'])
out = outdir / "{tsv_n}"
fields = ['transcript_id','n_positions','n_unique_back','n_multi_back','n_unmapped','uniqueness_factor']
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\\t')
    w.writeheader(); w.writerows(records)

fs = [r['uniqueness_factor'] for r in records]
n  = len(fs)
print(f"\\n  ── Distribution ─────────────────────────────────────")
print(f"  n         : {{n:,}}")
print(f"  mean      : {{sum(fs)/n:.4f}}")
print(f"  f=1.00    : {{sum(1 for f in fs if f==1.0):,}}")
print(f"  f>=0.90   : {{sum(1 for f in fs if f>=0.90):,}}")
print(f"  f 0.10-0.90: {{sum(1 for f in fs if 0.10<=f<1.0):,}}")
print(f"  f<0.10    : {{sum(1 for f in fs if f<0.10):,}}")
print(f"  f=0.00    : {{sum(1 for f in fs if f==0.0):,}}")
print(f"\\n  Saved: {{out}}")
PYEOF

echo ""
echo "════════════════════════════════════════════════════════"
echo "Done: $OUTDIR/{tsv_n}"
echo "════════════════════════════════════════════════════════"
'''

# ── Build and write each script ────────────────────────────────────────────────
def make_script(sid, ref, src, mode, L):
    odir_n = outdir_name(ref, src, mode, L)
    tsv_n  = tsv_filename(ref, src, mode, L)

    parts = [section_header(sid, ref, src, mode, L)]

    if src == 'exon_flank':
        parts.append(section_chr_map())

    parts.append(section_build_index(ref))

    if src == 'exon_flank':
        parts.append(section_simulate_exon_flank(mode))
    else:
        parts.append(section_simulate_mane_rna(mode))

    parts.append(section_star_align(ref, mode, odir_n))
    parts.append(section_compute_uniqueness(mode, tsv_n))

    return ''.join(parts)


if __name__ == '__main__':
    created = []
    for sid, ref, src, mode, L in NEW_SCRIPTS:
        content  = make_script(sid, ref, src, mode, L)
        out_path = HERE / f'mappability_correction_{sid}.sh'
        out_path.write_text(content)
        out_path.chmod(out_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        created.append((sid, ref, src, mode, L))
        print(f"  {sid:6s}  {REF_LABELS[ref]:<15}  {SRC_LABELS[src]:<12}  {mode}  {L:>3} bp  →  mappability_correction_{sid}.sh")

    print(f"\nCreated {len(created)} scripts.")
