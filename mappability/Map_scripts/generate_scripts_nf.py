#!/usr/bin/env python3
"""
Generate the remaining "genomic, nf reference" mappability scripts —
the _nf counterparts of v1_b/c/d, v3_a-d, v4_a-d, v6_a/b, and v8_a.

mappability_correction_v1_a_nf.sh (exon-flank, PE, 100 bp) was hand-written
first as the reference pattern this generator replicates. Only the genomic/
nf reference arm is covered — dedup_cds and MANE-transcriptome references are
unaffected by the reference swap and are intentionally left alone.

Each script is built from sections depending on:
  src      : exon_flank | mane_rna | mane_rna_proper (v8 only)
  mode     : PE | SE
  read_len : 75 | 100 | 150 | 200

mane_rna-sourced scripts need spliced transcript sequences for the new
reference. nf-core explicitly recommends NOT downloading a separate
transcriptome FASTA and instead deriving it from genome + annotation (see
Ref/recommendations_nf.txt), so Step 0c below splices MANE_Select-tagged
transcripts from Homo_sapiens.GRCh38.116.gtf + the primary-assembly FASTA into
a new MANE_RNA_nf.fna, built once and shared by all mane_rna-sourced scripts
here — mirroring how the STAR index is built once and reused.

The nf reference (genome, GTF, and everything derived from them) lives in
Ref/ alongside the original NCBI RefSeq files — not a separate nfRef/
directory — distinguished by filename/extension (Homo_sapiens.GRCh38.* vs.
GCF_*) rather than by location.

Writes only new mappability_correction_v*_nf.sh files; never overwrites the
existing (non-_nf) scripts.
"""

import stat
from pathlib import Path

HERE = Path(__file__).resolve().parent

# ── Scripts to generate: (script_id, src, mode, read_len) ─────────────────────
NEW_SCRIPTS_NF = [
    ('v1_b_nf', 'exon_flank', 'SE', 100),
    ('v1_c_nf', 'mane_rna',   'SE', 100),
    ('v1_d_nf', 'mane_rna',   'PE', 100),
    ('v3_a_nf', 'exon_flank', 'PE',  75),
    ('v3_b_nf', 'exon_flank', 'SE',  75),
    ('v3_c_nf', 'mane_rna',   'SE',  75),
    ('v3_d_nf', 'mane_rna',   'PE',  75),
    ('v4_a_nf', 'exon_flank', 'PE', 150),
    ('v4_b_nf', 'exon_flank', 'SE', 150),
    ('v4_c_nf', 'mane_rna',   'SE', 150),
    ('v4_d_nf', 'mane_rna',   'PE', 150),
    ('v6_a_nf', 'mane_rna',   'SE', 200),
    ('v6_b_nf', 'mane_rna',   'PE', 200),
]
# v8_a_nf (mane_rna_proper, PE, 150) is generated separately below —
# different fragment-based simulation, matching mappability_v8_HPC.sh's v8_a leg.

SRC_LABELS = {'exon_flank': 'exon-flank', 'mane_rna': 'MANE RNA', 'mane_rna_proper': 'MANE RNA (proper PE)'}


# ── Naming helpers — mirror the exact filenames the old genomic-reference
#    scripts already use (see results/transcript_uniqueness_factors_genomic_*),
#    with "_nf" appended before the extension / directory suffix. ────────────
def outdir_suffix(src, mode, L):
    # Matches the real historical filenames in results/transcript_uniqueness_factors_genomic_*:
    # every combination includes the L-tag, EXCEPT mane_rna+PE+100bp (the old v1_d.sh
    # hardcoded "transcript_uniqueness_factors_genomic_RNA_PE.tsv" with no L100bp suffix).
    if L == 100 and src == 'mane_rna' and mode == 'PE':
        L_tag = ''
    else:
        L_tag = f'_L{L}bp'
    if src == 'exon_flank':
        return f'genomic{"" if mode == "PE" else "_SE"}{L_tag}'
    else:
        return f'genomic_RNA_{mode}{L_tag}'

def outdir_name(src, mode, L, proper=False):
    suffix = outdir_suffix(src, mode, L)
    if proper:
        suffix += '_proper'
    return 'mappability_' + suffix + '_nf'

def tsv_filename(src, mode, L, proper=False):
    suffix = outdir_suffix(src, mode, L)
    if proper:
        suffix += '_proper'
    return 'transcript_uniqueness_factors_' + suffix + '_nf.tsv'


# ── Shared preamble (paths, common to every _nf genomic script) ──────────────
def section_header(sid, src, mode, L, proper=False):
    src_l = SRC_LABELS['mane_rna_proper' if proper else src]
    return f'''#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — {sid} — NF REFERENCE
#
# Reads     : {src_l}
# Mode      : {mode}
# Read len  : {L} bp
#
# nf reference: Ensembl GRCh38.116 primary assembly + matching GTF (see
# mappability_correction_v1_a_nf.sh header for the full rationale). Nothing
# here overwrites the existing (non-_nf) pipeline outputs, indices, or scripts.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
RESULTS="$SCRIPT_DIR/../results"

REF="$SCRIPT_DIR/../../Ref"
GENOME_FA_NF="$REF/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
GTF_NF="$REF/Homo_sapiens.GRCh38.116.gtf"
GENOME_IDX_NF="$REF/STAR_index_genome_nf"
MANE_RNA_NF="$REF/MANE_RNA_nf.fna"

OUTDIR="$RESULTS/{outdir_name(src, mode, L, proper)}"
THREADS=64
READ_LEN={L}
FLANK=50

mkdir -p "$OUTDIR/simreads" "$OUTDIR/bam" "$OUTDIR/logs"

# ── Load modules (standalone script — same set as run_all_mappability.sh /
#    mappability_v8_HPC.sh / run_3prime_HPC.sh) ───────────────────────────────
module purge
module load bio/STAR/2.7.11b-GCC-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0
module load lang/Python/3.11.5-GCCcore-13.2.0

# STAR BAM sorting needs many open file descriptors
ulimit -n 65536

echo "════════════════════════════════════════════════════════"
echo "Mappability — {sid} (nf reference, {src_l}, {mode}, {L} bp)"
echo "READ_LEN=$READ_LEN  THREADS=$THREADS"
echo "════════════════════════════════════════════════════════"
'''


def section_chr_map():
    """Step 0 — expand chr_name_map.tsv (exon-flank scripts only; QC/crosswalk artifact, not used internally).

    results/chr_name_map.tsv is NOT produced by any script — the old
    v1_a.sh writes its chr->NC map to its own per-run OUTDIR, never to a
    shared results/ location. A shared copy may or may not exist (someone
    may have copied one there by hand). Rather than depend on that, this
    step rebuilds the base chr->NC map itself, from the OLD NCBI genome
    FASTA headers, using the exact same header-parsing logic as v1_a.sh's
    own Step 0, if the shared file isn't already there — so this never
    fails just because nobody happened to copy chr_name_map.tsv to this
    deployment's results/ directory."""
    return r'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 0: Expand chr_name_map.tsv -> chr_name_map_nf.tsv (does not touch original)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0: Expanding chr_name_map.tsv with Ensembl-style names..."
if [ -f "$RESULTS/chr_name_map_nf.tsv" ]; then
    echo "Step 0: chr_name_map_nf.tsv already exists, skipping."
else
if [ ! -f "$RESULTS/chr_name_map.tsv" ]; then
    echo "  $RESULTS/chr_name_map.tsv not found — rebuilding from the old NCBI genome FASTA"
    echo "  headers (same logic as v1_a.sh's own Step 0), instead of depending on a"
    echo "  pre-existing shared copy that no script actually guarantees exists."
    OLD_GENOME_FA_FOR_CHR_MAP="$REF/GCF_000001405.40_GRCh38.p14_genomic.fna"
python3 - << PYEOF
import re
from pathlib import Path

genome_fa = Path("$OLD_GENOME_FA_FOR_CHR_MAP")
out_map   = Path("$RESULTS/chr_name_map.tsv")

ncbi_map = {}
with open(genome_fa) as fh:
    for line in fh:
        if not line.startswith('>'):
            continue
        line = line.strip()
        name = line[1:].split()[0]
        rest = line[1:]
        m = re.search(r'chromosome (\w+)[,\s]', rest)
        if m and name.startswith('NC_'):
            ncbi_map[f"chr{m.group(1)}"] = name
        if 'mitochondrion' in rest.lower() or 'mitochondrial' in rest.lower():
            ncbi_map['chrM'] = name
            ncbi_map['chrMT'] = name

with open(out_map, 'w') as fh:
    for k, v in sorted(ncbi_map.items()):
        fh.write(f"{k}\t{v}\n")
print(f"  Chromosomes mapped: {len(ncbi_map)}")
print(f"  Wrote: {out_map}")
PYEOF
fi
python3 - << PYEOF
from pathlib import Path

old_map_path = Path("$RESULTS/chr_name_map.tsv")
genome_fa_nf = Path("$GENOME_FA_NF")
out_map      = Path("$RESULTS/chr_name_map_nf.tsv")

def ensembl_name_for(chr_name):
    if chr_name in ("chrM", "chrMT"):
        return "MT"
    if chr_name.startswith("chr"):
        return chr_name[3:]
    return chr_name

rows = []
with open(old_map_path) as fh:
    for line in fh:
        parts = line.strip().split('\t')
        if len(parts) != 2:
            continue
        chr_name, ncbi_name = parts
        rows.append([chr_name, ncbi_name, ensembl_name_for(chr_name), 'canonical'])

canonical_ensembl_names = {r[2] for r in rows}
seen = set()
with open(genome_fa_nf) as fh:
    for line in fh:
        if line.startswith('>'):
            seen.add(line[1:].split()[0])

for name in sorted(seen - canonical_ensembl_names):
    rows.append(['NA', 'NA', name, 'unplaced_scaffold_needs_crosswalk'])

rows.sort(key=lambda r: (r[3] != 'canonical', r[0]))
with open(out_map, 'w') as fh:
    fh.write("chr_name\tncbi_name\tensembl_name\tnote\n")
    for r in rows:
        fh.write('\t'.join(r) + '\n')
print(f"  Wrote: {out_map}")
PYEOF
fi
echo "Step 0 complete."
'''


def section_build_index():
    """Step 0b — build STAR_index_genome_nf (shared; skipped if already exists)."""
    return r'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 0b: Build STAR index from nf genome (shared across all nf genomic scripts)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0b: Checking nf STAR index..."
mkdir -p "$GENOME_IDX_NF"
if [ ! -f "$GENOME_IDX_NF/SA" ]; then
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$GENOME_IDX_NF" \
        --genomeFastaFiles    "$GENOME_FA_NF" \
        --genomeSAindexNbases 14 \
        --genomeChrBinNbits   14 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR/logs/star_index_nf.log"
    echo "Step 0b complete."
else
    echo "Step 0b: nf index already exists, skipping."
fi
'''


def section_build_mane_rna():
    """Step 0c — splice MANE_Select transcripts from GTF_NF + genome into MANE_RNA_nf.fna
    (shared across all mane_rna-sourced nf scripts; skipped if already built)."""
    return r'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 0c: Build MANE_RNA_nf.fna — spliced MANE_Select transcript sequences,
# derived from the nf genome + GTF (nf-core recommends deriving this rather
# than downloading a separate transcriptome FASTA). Shared; skipped if exists.
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0c: Checking MANE_RNA_nf.fna..."
if [ -f "$MANE_RNA_NF" ]; then
    echo "Step 0c: MANE_RNA_nf.fna already exists, skipping."
else
    echo "Step 0c: Splicing MANE_Select transcripts from nf GTF..."
    if [ ! -f "$GENOME_FA_NF.fai" ]; then
        samtools faidx "$GENOME_FA_NF"
    fi
python3 - << PYEOF
import subprocess
from pathlib import Path
from collections import defaultdict

gtf_path  = Path("$GTF_NF")
genome_fa = Path("$GENOME_FA_NF")
out_fa    = Path("$MANE_RNA_NF")

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def revcomp(s): return s.translate(RC)[::-1]

def parse_gtf_attrs(col9):
    attrs = defaultdict(list)
    for field in col9.strip().split(';'):
        field = field.strip()
        if not field:
            continue
        key, _, val = field.partition(' ')
        attrs[key].append(val.strip('"'))
    return attrs

mane_select_tx = set()
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'transcript':
            continue
        a = parse_gtf_attrs(f[8])
        if 'MANE_Select' not in a.get('tag', []):
            continue
        tid  = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None:
            continue
        mane_select_tx.add(f"{tid}.{tver}" if tver else tid)
print(f"  MANE_Select transcripts: {len(mane_select_tx):,}", flush=True)

tx_exons = defaultdict(list)
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        a = parse_gtf_attrs(f[8])
        tid  = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None:
            continue
        versioned_tid = f"{tid}.{tver}" if tver else tid
        if versioned_tid not in mane_select_tx:
            continue
        tx_exons[versioned_tid].append((f[0], int(f[3]) - 1, int(f[4]), f[6]))

for tid in tx_exons:
    tx_exons[tid].sort(key=lambda e: e[1])

n_written = n_skip = 0
with open(out_fa, 'w') as out:
    for tid, exons in tx_exons.items():
        chrom  = exons[0][0]
        strand = exons[0][3]
        seq_parts = []
        ok = True
        for (c, start, end, _s) in exons:
            region = f"{c}:{start + 1}-{end}"
            try:
                res = subprocess.run(['samtools', 'faidx', str(genome_fa), region],
                                     capture_output=True, text=True, check=True)
                seq_parts.append(''.join(res.stdout.split('\n')[1:]).upper())
            except subprocess.CalledProcessError:
                ok = False
                break
        if not ok:
            n_skip += 1
            continue
        seq = ''.join(seq_parts)
        if strand == '-':
            seq = revcomp(seq)
        out.write(f">{tid}\n{seq}\n")
        n_written += 1

print(f"  Transcripts written: {n_written:,}  (skipped: {n_skip})")
print(f"  Saved: {out_fa}", flush=True)
PYEOF
fi
echo "Step 0c complete."
'''


def section_simulate_exon_flank(mode):
    """Step 1 — exon-flank read simulation from GTF_NF (MANE_Select filter) + nf genome."""
    if mode == 'PE':
        open_files   = "with open(r1_path, 'w') as r1, open(r2_path, 'w') as r2:"
        r2_decl      = "r2_path = outdir / 'simreads/sim_R2.fastq'"
        write_reads  = (
            "                r1.write(f'@{name}\\n{fwd}\\n+\\n{qual}\\n')\n"
            "                r2.write(f'@{name}\\n{revcomp(fwd)}\\n+\\n{qual}\\n')"
        )
        revcomp_func = "RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')\ndef revcomp(s): return s.translate(RC)[::-1]\n"
    else:
        open_files   = "with open(r1_path, 'w') as r1:"
        r2_decl      = ''
        write_reads  = "                r1.write(f'@{name}\\n{fwd}\\n+\\n{qual}\\n')"
        revcomp_func = ''

    return f'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Parse nf GTF (MANE_Select filter), fetch exon ± FLANK, simulate reads
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Simulating exon-flank reads (READ_LEN=$READ_LEN, FLANK=$FLANK, {mode})..."

python3 - << PYEOF
import subprocess
from pathlib import Path
from collections import defaultdict

gtf_path  = Path("$GTF_NF")
genome_fa = Path("$GENOME_FA_NF")
outdir    = Path("$OUTDIR")
read_len  = int("$READ_LEN")
flank     = int("$FLANK")

{revcomp_func}
def parse_gtf_attrs(col9):
    attrs = defaultdict(list)
    for field in col9.strip().split(';'):
        field = field.strip()
        if not field: continue
        key, _, val = field.partition(' ')
        attrs[key].append(val.strip('"'))
    return attrs

mane_select_tx = set()
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip('\\n').split('\\t')
        if len(f) < 9 or f[2] != 'transcript': continue
        a = parse_gtf_attrs(f[8])
        if 'MANE_Select' not in a.get('tag', []): continue
        tid = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None: continue
        mane_select_tx.add(f"{{tid}}.{{tver}}" if tver else tid)
print(f"  MANE_Select transcripts: {{len(mane_select_tx):,}}", flush=True)

tx_exons = defaultdict(list)
n_exons = 0
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip('\\n').split('\\t')
        if len(f) < 9 or f[2] != 'exon': continue
        a = parse_gtf_attrs(f[8])
        tid = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None: continue
        versioned_tid = f"{{tid}}.{{tver}}" if tver else tid
        if versioned_tid not in mane_select_tx: continue
        chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
        tx_exons[versioned_tid].append((chrom, start, end, strand))
        n_exons += 1
print(f"  Exon features parsed: {{n_exons:,}}", flush=True)

fai = Path(str(genome_fa) + '.fai')
if not fai.exists():
    subprocess.run(['samtools', 'faidx', str(genome_fa)], check=True)

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
    """Step 1 — read simulation from spliced MANE_RNA_nf.fna."""
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
# STEP 1: Simulate reads from spliced MANE_RNA_nf.fna
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Simulating MANE RNA reads (READ_LEN=$READ_LEN, {mode})..."

python3 - << PYEOF
from pathlib import Path

mane_rna = Path("$MANE_RNA_NF")
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


def section_simulate_mane_rna_proper(read_len, frag_len):
    """Step 1 — v8-style proper PE fragment simulation from spliced MANE_RNA_nf.fna."""
    return f'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Simulate proper PE reads from spliced MANE_RNA_nf.fna
# (fragment = {frag_len} bp; R1 = first {read_len} bp fwd; R2 = last {read_len} bp revcomp)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Simulating proper PE MANE RNA reads (READ_LEN=$READ_LEN, FRAG_LEN={frag_len})..."

python3 - << PYEOF
from pathlib import Path

mane_rna = Path("$MANE_RNA_NF")
outdir   = Path("$OUTDIR")
read_len = int("$READ_LEN")
frag_len = {frag_len}

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def revcomp(s): return s.translate(RC)[::-1]

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

qual = 'I' * read_len
n_reads = n_skip = n_seq = 0

with open(outdir / 'simreads/sim_R1.fastq', 'w') as r1, \\
     open(outdir / 'simreads/sim_R2.fastq', 'w') as r2:
    for tid, seq in parse_fasta(mane_rna):
        L = len(seq)
        if L < frag_len: n_skip += 1; continue
        n_seq += 1
        for pos in range(L - frag_len + 1):
            fwd = seq[pos            : pos + read_len]
            rev = seq[pos + read_len : pos + frag_len]
            if fwd.count('N') > read_len // 5: continue
            if rev.count('N') > read_len // 5: continue
            r1.write(f'@{{tid}}|{{pos}}\\n{{fwd}}\\n+\\n{{qual}}\\n')
            r2.write(f'@{{tid}}|{{pos}}\\n{{revcomp(rev)}}\\n+\\n{{qual}}\\n')
            n_reads += 1

print(f"  Transcripts: {{n_seq:,}}  (too short for {{frag_len}} bp: {{n_skip}})")
print(f"  Read pairs: {{n_reads:,}}", flush=True)
PYEOF
echo "Step 1 complete."
'''


def section_star_align(mode):
    """Step 2 — STAR alignment against the nf genome index."""
    if mode == 'PE':
        reads_in = (
            '    --readFilesIn           "$OUTDIR/simreads/sim_R1.fastq" \\\n'
            '                            "$OUTDIR/simreads/sim_R2.fastq" \\'
        )
    else:
        reads_in = '    --readFilesIn           "$OUTDIR/simreads/sim_R1.fastq" \\'

    return f'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: STAR alignment against nf genome index — {mode}
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 2: STAR alignment ({mode}, nf genome, $THREADS threads)..."

STAR \\
    --runThreadN            $THREADS \\
    --genomeDir             "$GENOME_IDX_NF" \\
{reads_in}
    --outSAMtype            BAM SortedByCoordinate \\
    --outSAMattributes      NH HI AS NM \\
    --outSAMmultNmax        1 \\
    --outFilterMultimapNmax 20 \\
    --outBAMsortingThreadN  $THREADS \\
    --outBAMsortingBinsN    20 \\
    --limitBAMsortRAM       160000000000 \\
    --outFileNamePrefix     "$OUTDIR/bam/sim_" \\
    2>&1 | tee "$OUTDIR/logs/star_sim_nf.log"

samtools index -@ $THREADS "$OUTDIR/bam/sim_Aligned.sortedByCoord.out.bam"
echo "Step 2 complete."
'''


def section_compute_uniqueness(mode, tsv_n):
    """Step 3 — compute per-transcript uniqueness factor."""
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
import pysam, csv
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

fs = [r['uniqueness_factor'] for r in records]; n = len(fs)
print(f"\\n  ── Distribution ─────────────────────────────────────")
print(f"  n         : {{n:,}}")
if n:
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
echo "Compare against: results/{tsv_n.replace('_nf.tsv', '.tsv')}"
echo "════════════════════════════════════════════════════════"
'''


# ── Assemble a full script ─────────────────────────────────────────────────────
def make_script(sid, src, mode, L):
    parts = [section_header(sid, src, mode, L)]
    if src == 'exon_flank':
        parts.append(section_chr_map())
    parts.append(section_build_index())
    if src == 'mane_rna':
        parts.append(section_build_mane_rna())
    if src == 'exon_flank':
        parts.append(section_simulate_exon_flank(mode))
    else:
        parts.append(section_simulate_mane_rna(mode))
    parts.append(section_star_align(mode))
    parts.append(section_compute_uniqueness(mode, tsv_filename(src, mode, L)))
    return ''.join(parts)


def make_script_v8a_nf():
    sid, mode, L, frag = 'v8_a_nf', 'PE', 150, 300
    parts = [section_header(sid, 'mane_rna', mode, L, proper=True)]
    parts.append(section_build_index())
    parts.append(section_build_mane_rna())
    parts.append(section_simulate_mane_rna_proper(L, frag))
    parts.append(section_star_align(mode))
    parts.append(section_compute_uniqueness(mode, tsv_filename('mane_rna', mode, L, proper=True)))
    return ''.join(parts)


# ═══════════════════════════════════════════════════════════════════════════
# v7_nf / v8_b_nf / v8_c_nf — dedup-CDS reference (genome-wide, protein-coding)
#
# Unlike the MANE_Select-only scope used everywhere above, the old v2/v7/v8_b
# dedup-CDS reference covers EVERY protein-coding transcript in the genome
# (93,088 sequences after seqkit dedup) — deliberately not MANE-restricted, so
# it is a fair like-for-like comparison against the old v7/v8_b results.
# ═══════════════════════════════════════════════════════════════════════════

NEW_SCRIPTS_V7_NF = [
    ('v7_a_nf', 'SE', 100, 'genome'),
    ('v7_b_nf', 'PE', 100, 'genome'),
    ('v7_c_nf', 'SE',  75, 'genome'),
    ('v7_d_nf', 'PE',  75, 'genome'),
    ('v7_e_nf', 'SE', 150, 'genome'),
    ('v7_f_nf', 'PE', 150, 'genome'),
    ('v7_g_nf', 'SE', 100, 'selfmap'),
    ('v7_h_nf', 'PE', 100, 'selfmap'),
    ('v7_i_nf', 'SE', 200, 'genome'),
    ('v7_j_nf', 'PE', 200, 'genome'),
]


def v7_outdir_name(mode, L, target):
    if target == 'selfmap':
        return f'mappability_dedup_cds_reads_selfmap_{mode}_nf'
    L_tag = '' if L == 100 else f'_L{L}bp'
    return f'mappability_dedup_cds_reads_{mode}{L_tag}_nf'

def v7_tsv_name(mode, L, target):
    if target == 'selfmap':
        return f'cds_uniqueness_factors_selfmap_{mode}_nf.tsv'
    L_tag = '' if L == 100 else f'_L{L}bp'
    return f'cds_uniqueness_factors_genome_{mode}{L_tag}_nf.tsv'


def section_header_v7(sid, mode, L, target):
    label = 'full genome' if target == 'genome' else 'dedup CDS self-map'
    return f'''#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — {sid} — NF REFERENCE
#
# Read source : deduplicated CDS, genome-wide protein-coding (nf reference)
# Reference   : {label}
# Mode        : {mode}
# Read len    : {L} bp
#
# nf counterpart of the old v7 family. See mappability_correction_v1_a_nf.sh
# for the general nf-reference rationale, and Step 0d/0e below for how the
# genome-wide dedup-CDS reference is built for the new genome + GTF.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
RESULTS="$SCRIPT_DIR/../results"

REF="$SCRIPT_DIR/../../Ref"
GENOME_FA_NF="$REF/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
GTF_NF="$REF/Homo_sapiens.GRCh38.116.gtf"
GENOME_IDX_NF="$REF/STAR_index_genome_nf"
DEDUP_CDS_NF_RAW="$REF/Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_raw.fna"
DEDUP_CDS_NF="$REF/Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna"
CDS_IDX_NF="$REF/STAR_index_dedup_cds_nf"

OUTDIR="$RESULTS/{v7_outdir_name(mode, L, target)}"
THREADS=64
READ_LEN={L}

mkdir -p "$OUTDIR/simreads" "$OUTDIR/bam" "$OUTDIR/logs"

module purge
module load bio/STAR/2.7.11b-GCC-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0
module load lang/Python/3.11.5-GCCcore-13.2.0

ulimit -n 65536

echo "════════════════════════════════════════════════════════"
echo "Mappability — {sid} (nf reference, dedup CDS reads, {label}, {mode}, {L} bp)"
echo "READ_LEN=$READ_LEN  THREADS=$THREADS"
echo "════════════════════════════════════════════════════════"
'''


def section_build_dedup_cds_nf():
    """Step 0d — extract genome-wide protein-coding CDS from GTF_NF + genome,
    then dedup with seqkit. Shared across all v7_nf/v8_b_nf scripts; skipped
    if Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna already exists. Uses pysam.FastaFile (not per-exon
    samtools subprocess calls) since this covers far more transcripts than
    the MANE_Select-only extraction elsewhere."""
    return r'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 0d: Build Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna — ALL protein-coding CDS, genome-wide,
# deduplicated with seqkit (shared; skipped if already built)
#
# Split into two independently-resumable halves because seqkit has no
# confirmed module on this cluster: extraction (pure Python + pysam, no
# extra dependency) always runs here; dedup only runs here if seqkit is
# already on PATH, otherwise it stops with instructions to scp the raw
# FASTA down, run 'seqkit rmdup -s' locally, and scp Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna back
# to the same path here — the next run then just skips straight past this.
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0d: Checking Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna..."
if [ -f "$DEDUP_CDS_NF" ]; then
    echo "Step 0d: Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna already exists, skipping."
    SKIP_STEP_0D=1
else
    SKIP_STEP_0D=0
fi

if [ "$SKIP_STEP_0D" = "0" ] && [ -f "$DEDUP_CDS_NF_RAW" ]; then
    echo "Step 0d: Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_raw.fna already extracted, skipping extraction."
elif [ "$SKIP_STEP_0D" = "0" ]; then
    echo "Step 0d: Extracting genome-wide protein-coding CDS from nf GTF..."
python3 - << PYEOF
import pysam
from pathlib import Path
from collections import defaultdict

gtf_path  = Path("$GTF_NF")
genome_fa = Path("$GENOME_FA_NF")
raw_fa    = Path("$DEDUP_CDS_NF_RAW")

def parse_gtf_attrs(col9):
    attrs = defaultdict(list)
    for field in col9.strip().split(';'):
        field = field.strip()
        if not field:
            continue
        key, _, val = field.partition(' ')
        attrs[key].append(val.strip('"'))
    return attrs

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def revcomp(s): return s.translate(RC)[::-1]

print("  Scanning GTF for protein-coding CDS features...", flush=True)
tx_cds = defaultdict(list)
n_cds_lines = 0
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        a = parse_gtf_attrs(f[8])
        biotype = (a.get('transcript_biotype') or a.get('gene_biotype') or [''])[0]
        if biotype != 'protein_coding':
            continue
        tid = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None:
            continue
        versioned_tid = f"{tid}.{tver}" if tver else tid
        chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
        tx_cds[versioned_tid].append((chrom, start, end, strand))
        n_cds_lines += 1

for tid in tx_cds:
    tx_cds[tid].sort(key=lambda e: e[1])

print(f"  CDS feature lines: {n_cds_lines:,}")
print(f"  Protein-coding transcripts with CDS: {len(tx_cds):,}", flush=True)

print("  Fetching + splicing CDS sequences...", flush=True)
fa = pysam.FastaFile(str(genome_fa))
n_written = n_skip = 0
with open(raw_fa, 'w') as out:
    for tid, segs in tx_cds.items():
        strand = segs[0][3]
        try:
            seq = ''.join(fa.fetch(c, s, e) for (c, s, e, _st) in segs).upper()
        except (KeyError, ValueError):
            n_skip += 1
            continue
        if not seq:
            n_skip += 1
            continue
        if strand == '-':
            seq = revcomp(seq)
        out.write(f">{tid}\n{seq}\n")
        n_written += 1
fa.close()

print(f"  CDS sequences written (pre-dedup): {n_written:,}  (skipped: {n_skip})")
print(f"  Saved raw: {raw_fa}", flush=True)
PYEOF
fi

if [ "$SKIP_STEP_0D" = "0" ]; then
    if command -v seqkit &> /dev/null; then
        echo "  Deduplicating exact-sequence duplicates with seqkit rmdup -s..."
        seqkit rmdup -s -i "$DEDUP_CDS_NF_RAW" -o "$DEDUP_CDS_NF" --threads $THREADS
        echo "Step 0d complete."
    else
        echo ""
        echo "  seqkit not found on PATH on this cluster."
        echo "  $DEDUP_CDS_NF_RAW is ready (genome-wide protein-coding CDS, not yet deduplicated)."
        echo "  To finish Step 0d:"
        echo "    1. scp $DEDUP_CDS_NF_RAW to a machine with seqkit installed"
        echo "    2. run: seqkit rmdup -s -i Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_raw.fna -o Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna"
        echo "    3. scp Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna back to: $DEDUP_CDS_NF"
        echo "    4. re-run this script — it will detect the file and skip straight past Step 0d"
        exit 1
    fi
fi
'''


def section_build_cds_index_nf():
    """Step 0e — build STAR_index_dedup_cds_nf (shared; skipped if exists)."""
    return r'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 0e: Build STAR index (dedup CDS, nf) — shared; skipped if already exists
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0e: Checking STAR dedup-CDS-nf index..."
mkdir -p "$CDS_IDX_NF"
if [ ! -f "$CDS_IDX_NF/SA" ]; then
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$CDS_IDX_NF" \
        --genomeFastaFiles    "$DEDUP_CDS_NF" \
        --genomeSAindexNbases 12 \
        --genomeChrBinNbits   11 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR/logs/star_index_cds_nf.log"
    echo "Step 0e complete."
else
    echo "Step 0e: index already exists, skipping."
fi
'''


def section_simulate_cds_reads_nf(mode):
    """Step 1 — simulate reads from Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna (SE or fake-PE)."""
    if mode == 'PE':
        open_files   = "with open(r1_path, 'w') as r1, open(r2_path, 'w') as r2:"
        r2_decl      = "r2_path = outdir / 'simreads/sim_R2.fastq'"
        write_reads  = (
            "            r1.write(f'@{seq_id}|{pos}\\n{fwd}\\n+\\n{qual}\\n')\n"
            "            r2.write(f'@{seq_id}|{pos}\\n{revcomp(fwd)}\\n+\\n{qual}\\n')"
        )
        revcomp_func = "RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')\ndef revcomp(s): return s.translate(RC)[::-1]\n"
    else:
        open_files   = "with open(r1_path, 'w') as r1:"
        r2_decl      = ''
        write_reads  = "            r1.write(f'@{seq_id}|{pos}\\n{fwd}\\n+\\n{qual}\\n')"
        revcomp_func = ''

    return f'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Simulate reads from Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna ({mode})
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Simulating reads from Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna (READ_LEN=$READ_LEN, {mode})..."

if [ -f "$OUTDIR/simreads/sim_R1.fastq" ]; then
    echo "  Reads exist, skipping simulation."
else

python3 - << PYEOF
from pathlib import Path

cds_fa   = Path("$DEDUP_CDS_NF")
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
n_reads = n_seq = n_skip = 0

{open_files}
    for seq_id, seq in parse_fasta(cds_fa):
        L = len(seq)
        if L < read_len: n_skip += 1; continue
        n_seq += 1
        for pos in range(L - read_len + 1):
            fwd = seq[pos : pos + read_len]
            if fwd.count('N') > read_len // 5: continue
{write_reads}
            n_reads += 1

print(f"  CDS sequences: {{n_seq:,}}  (too short: {{n_skip}})")
print(f"  Total reads simulated: {{n_reads:,}}", flush=True)
PYEOF

fi
echo "Step 1 complete."
'''


def section_star_align_v7(mode, target):
    """Step 2 — STAR alignment against genome_nf or dedup_cds_nf (self-map)."""
    idx = '"$GENOME_IDX_NF"' if target == 'genome' else '"$CDS_IDX_NF"'
    label = 'full genome' if target == 'genome' else 'dedup CDS self-map'
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
    --outFilterMultimapNmax 20 \\
    --outBAMsortingThreadN  $THREADS \\
    --outBAMsortingBinsN    20 \\
    --limitBAMsortRAM       160000000000 \\
    --outFileNamePrefix     "$OUTDIR/bam/sim_" \\
    2>&1 | tee "$OUTDIR/logs/star_sim.log"

samtools index -@ $THREADS "$OUTDIR/bam/sim_Aligned.sortedByCoord.out.bam"
echo "Step 2 complete."
'''


def section_compute_uniqueness_cds(mode, tsv_n):
    """Step 3 — compute per-CDS uniqueness factor. CDS IDs here are plain
    versioned transcript IDs with no embedded '|', so unlike the old v7
    scripts (which had to be patched post-hoc — see mappability_v7_fix_uf.sh
    — because NCBI's 'lcl|...' IDs contain '|' themselves), a plain
    split('|')[0] on the read name correctly recovers the full CDS ID."""
    if mode == 'PE':
        bam_filter = 'if read.is_supplementary or read.is_secondary or read.is_read2: continue'
    else:
        bam_filter = 'if read.is_supplementary or read.is_secondary or read.is_unmapped: continue'

    return f'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Compute per-CDS uniqueness factor
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 3: Computing CDS uniqueness factors..."

python3 - << PYEOF
import pysam, csv
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
print(f"  CDS seqs: {{len(sim_counts):,}}  Reads: {{sum(sim_counts.values()):,}}", flush=True)

unique_back = defaultdict(int)
total_back  = defaultdict(int)
bam = pysam.AlignmentFile(str(bam_path), 'rb')
for read in bam.fetch():
    {bam_filter}
    cid = read.query_name.split('|')[0]
    total_back[cid] += 1
    if read.mapping_quality == 255: unique_back[cid] += 1
bam.close()

records = []
for cid, n_sim in sim_counts.items():
    n_u = unique_back.get(cid, 0); n_t = total_back.get(cid, 0)
    records.append({{
        'cds_id':            cid,
        'n_positions':       n_sim,
        'n_unique_back':     n_u,
        'n_multi_back':      n_t - n_u,
        'n_unmapped':        n_sim - n_t,
        'uniqueness_factor': round(n_u / n_sim, 6) if n_sim > 0 else 0.0,
    }})

records.sort(key=lambda r: r['uniqueness_factor'])
out = outdir / "{tsv_n}"
fields = ['cds_id','n_positions','n_unique_back','n_multi_back','n_unmapped','uniqueness_factor']
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\\t')
    w.writeheader(); w.writerows(records)

fs = [r['uniqueness_factor'] for r in records]; n = len(fs)
print(f"\\n  ── Distribution ─────────────────────────────────────")
print(f"  n          : {{n:,}}")
if n:
    print(f"  mean       : {{sum(fs)/n:.4f}}")
    print(f"  f=1.00     : {{sum(1 for f in fs if f==1.0):,}}")
    print(f"  f>=0.90    : {{sum(1 for f in fs if f>=0.90):,}}")
    print(f"  f<0.10     : {{sum(1 for f in fs if f<0.10):,}}")
    print(f"  f=0.00     : {{sum(1 for f in fs if f==0.0):,}}")
print(f"\\n  Saved: {{out}}")
PYEOF

echo ""
echo "════════════════════════════════════════════════════════"
echo "Done: $OUTDIR/{tsv_n}"
echo "Compare against: results/{tsv_n.replace('_nf.tsv', '.tsv')}"
echo "════════════════════════════════════════════════════════"
'''


def make_script_v7(sid, mode, L, target):
    parts = [section_header_v7(sid, mode, L, target)]
    parts.append(section_build_dedup_cds_nf())
    parts.append(section_build_cds_index_nf())
    if target == 'genome':
        parts.append(section_build_index())   # STAR_index_genome_nf (shared w/ v1-v6_nf)
    parts.append(section_simulate_cds_reads_nf(mode))
    parts.append(section_star_align_v7(mode, target))
    parts.append(section_compute_uniqueness_cds(mode, v7_tsv_name(mode, L, target)))
    return ''.join(parts)


# ── v8_b_nf / v8_c_nf — reuse v8_a_nf's simulated reads ───────────────────────

def make_script_v8bc_nf(which):
    """which: 'b' (dedup CDS) or 'c' (MANE RNA self-map)."""
    if which == 'b':
        sid, ref_label, idx_var, idx_build, tsv = (
            'v8_b_nf', 'dedup CDS (nf)', 'CDS_IDX_NF',
            section_build_cds_index_nf(),
            'transcript_uniqueness_factors_dedup_cds_RNA_PE_L150bp_proper_nf.tsv',
        )
        extra_index_deps = section_build_dedup_cds_nf() + idx_build
        outdir = 'mappability_dedup_cds_RNA_PE_L150bp_proper_nf'
        idx_path = '"$CDS_IDX_NF"'
    else:
        sid, ref_label, tsv = (
            'v8_c_nf', 'MANE RNA (nf, self-map)',
            'transcript_uniqueness_factors_MANE_RNA_PE_L150bp_proper_nf.tsv',
        )
        outdir = 'mappability_MANE_RNA_PE_L150bp_proper_nf'
        idx_path = '"$MANE_IDX_NF"'
        extra_index_deps = r'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 0: Build STAR index (MANE RNA, nf) — shared; skipped if already exists
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0: Checking STAR MANE-RNA-nf index..."
mkdir -p "$MANE_IDX_NF"
if [ ! -f "$MANE_IDX_NF/SA" ]; then
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$MANE_IDX_NF" \
        --genomeFastaFiles    "$MANE_RNA_NF" \
        --genomeSAindexNbases 12 \
        --genomeChrBinNbits   11 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR/logs/star_index_mane_nf.log"
else
    echo "  Index exists, skipping."
fi
'''

    header = f'''#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — {sid} — NF REFERENCE
#
# Reuses reads simulated by mappability_correction_v8_a_nf.sh (MANE RNA,
# proper PE, 150 bp) — run v8_a_nf first. Aligns against {ref_label}
# instead of the nf genome, to isolate the competitor-set effect under
# realistic PE geometry (same comparison as the old v8_b/v8_c).
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
RESULTS="$SCRIPT_DIR/../results"

REF="$SCRIPT_DIR/../../Ref"
GENOME_FA_NF="$REF/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
GTF_NF="$REF/Homo_sapiens.GRCh38.116.gtf"
MANE_RNA_NF="$REF/MANE_RNA_nf.fna"
MANE_IDX_NF="$REF/STAR_index_MANE_nf"
DEDUP_CDS_NF_RAW="$REF/Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_raw.fna"
DEDUP_CDS_NF="$REF/Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna"
CDS_IDX_NF="$REF/STAR_index_dedup_cds_nf"

SIM_DIR="$RESULTS/mappability_genomic_RNA_PE_L150bp_proper_nf/simreads"
OUTDIR="$RESULTS/{outdir}"
THREADS=64

mkdir -p "$OUTDIR/bam" "$OUTDIR/logs"

module purge
module load bio/STAR/2.7.11b-GCC-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0
module load lang/Python/3.11.5-GCCcore-13.2.0

ulimit -n 65536

echo "════════════════════════════════════════════════════════"
echo "Mappability — {sid} (MANE RNA reads, proper PE, 150 bp, {ref_label})"
echo "Reads from: $SIM_DIR"
echo "════════════════════════════════════════════════════════"
'''

    body = f'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Check reads from v8_a_nf
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 1: Checking for reads from v8_a_nf..."
if [ ! -f "$SIM_DIR/sim_R1.fastq" ]; then
    echo "ERROR: reads not found at $SIM_DIR — run v8_a_nf first." >&2
    exit 1
fi
echo "  Reads found, proceeding."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: STAR alignment (PE, {ref_label})
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 2: STAR alignment (PE, {ref_label}, $THREADS threads)..."

STAR \\
    --runThreadN            $THREADS \\
    --genomeDir             {idx_path} \\
    --readFilesIn           "$SIM_DIR/sim_R1.fastq" \\
                            "$SIM_DIR/sim_R2.fastq" \\
    --outSAMtype            BAM SortedByCoordinate \\
    --outSAMattributes      NH HI AS NM \\
    --outSAMmultNmax        1 \\
    --outFilterMultimapNmax 20 \\
    --alignIntronMax        1 \\
    --outBAMsortingThreadN  $THREADS \\
    --outBAMsortingBinsN    20 \\
    --limitBAMsortRAM       160000000000 \\
    --outFileNamePrefix     "$OUTDIR/bam/sim_" \\
    2>&1 | tee "$OUTDIR/logs/star_sim.log"

samtools index -@ $THREADS "$OUTDIR/bam/sim_Aligned.sortedByCoord.out.bam"
echo "Step 2 complete."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Compute transcript uniqueness factors
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 3: Computing transcript uniqueness factors..."

python3 - << PYEOF
import pysam, csv
from collections import defaultdict
from pathlib import Path

simdir   = Path("$SIM_DIR")
outdir   = Path("$OUTDIR")
bam_path = outdir / "bam/sim_Aligned.sortedByCoord.out.bam"
r1_path  = simdir / "sim_R1.fastq"

sim_counts = defaultdict(int)
with open(r1_path) as fh:
    for line in fh:
        if line.startswith('@'):
            sim_counts[line[1:].split('|')[0]] += 1

unique_back = defaultdict(int)
total_back  = defaultdict(int)
bam = pysam.AlignmentFile(str(bam_path), 'rb')
for read in bam.fetch():
    if read.is_supplementary or read.is_secondary or read.is_read2: continue
    tid = read.query_name.split('|')[0]
    total_back[tid] += 1
    if read.mapping_quality == 255: unique_back[tid] += 1
bam.close()

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
out = outdir / "{tsv}"
fields = ['transcript_id','n_positions','n_unique_back','n_multi_back','n_unmapped','uniqueness_factor']
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\\t')
    w.writeheader(); w.writerows(records)

fs = [r['uniqueness_factor'] for r in records]; n = len(fs)
print(f"  n={{n:,}}  mean={{sum(fs)/n:.4f}}  UF=1.0: {{sum(1 for f in fs if f==1.0):,}}  UF=0.0: {{sum(1 for f in fs if f==0.0):,}}")
print(f"  Saved: {{out}}")
PYEOF

echo "════════════════════════════════════════════════════════"
echo "Done: $OUTDIR/{tsv}"
echo "════════════════════════════════════════════════════════"
'''
    return header + extra_index_deps + body


# ═══════════════════════════════════════════════════════════════════════════
# v2_nf — dedup-CDS reference, MANE_Select-restricted (NOT the genome-wide
# reference built for v7_nf/v8_b_nf above)
#
# The old v2 scripts align MANE-sourced reads (exon-flank / MANE RNA) against
# the SAME genome-wide dedup_cds.fna used by v7 (confirmed: v2_a.sh/v2_c.sh
# both point at STAR_index_dedup_cds). For the nf counterpart, a separate,
# smaller MANE_Select-only CDS reference is used instead — a third reference
# variant, distinct from both old v2 and genome-wide Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna.
# ═══════════════════════════════════════════════════════════════════════════

# (sid, src, mode, L, align_intron_max1, outdir_name, tsv_name) — exact names
# mirror the old v2 scripts' idiosyncratic naming (e.g. v2_c's bare
# "dedup_cds" outdir, v2_e's missing L100bp tag), with "_nf" appended.
OLD_V2 = {
    'v2_a': ('exon_flank', 'PE', 100, False, 'mappability_dedup_cds_PE_flank', 'transcript_uniqueness_factors_dedup_cds_PE_flank_L100bp.tsv'),
    'v2_b': ('exon_flank', 'SE', 100, False, 'mappability_dedup_cds_SE_flank', 'transcript_uniqueness_factors_dedup_cds_SE_flank_L100bp.tsv'),
    'v2_c': ('mane_rna',   'SE', 100, False, 'mappability_dedup_cds',          'transcript_uniqueness_factors_dedup_cds_L100bp.tsv'),
    'v2_d': ('mane_rna',   'SE', 100, True,  'mappability_dedup_cds_alignIntronMax1', 'transcript_uniqueness_factors_dedup_cds_alignIntronMax1_L100bp.tsv'),
    'v2_e': ('mane_rna',   'PE', 100, False, 'mappability_dedup_cds_RNA_PE',   'transcript_uniqueness_factors_dedup_cds_RNA_PE.tsv'),
    'v2_f': ('exon_flank', 'PE',  75, False, 'mappability_dedup_cds_PE_flank_L75bp', 'transcript_uniqueness_factors_dedup_cds_PE_flank_L75bp.tsv'),
    'v2_g': ('exon_flank', 'SE',  75, False, 'mappability_dedup_cds_SE_flank_L75bp', 'transcript_uniqueness_factors_dedup_cds_SE_flank_L75bp.tsv'),
    'v2_h': ('exon_flank', 'PE', 150, False, 'mappability_dedup_cds_PE_flank_L150bp', 'transcript_uniqueness_factors_dedup_cds_PE_flank_L150bp.tsv'),
    'v2_i': ('exon_flank', 'SE', 150, False, 'mappability_dedup_cds_SE_flank_L150bp', 'transcript_uniqueness_factors_dedup_cds_SE_flank_L150bp.tsv'),
    'v2_j': ('mane_rna',   'SE',  75, False, 'mappability_dedup_cds_RNA_SE_L75bp', 'transcript_uniqueness_factors_dedup_cds_RNA_SE_L75bp.tsv'),
    'v2_k': ('mane_rna',   'PE',  75, False, 'mappability_dedup_cds_RNA_PE_L75bp', 'transcript_uniqueness_factors_dedup_cds_RNA_PE_L75bp.tsv'),
    'v2_l': ('mane_rna',   'SE', 150, False, 'mappability_dedup_cds_RNA_SE_L150bp', 'transcript_uniqueness_factors_dedup_cds_RNA_SE_L150bp.tsv'),
    'v2_m': ('mane_rna',   'PE', 150, False, 'mappability_dedup_cds_RNA_PE_L150bp', 'transcript_uniqueness_factors_dedup_cds_RNA_PE_L150bp.tsv'),
    'v2_n': ('mane_rna',   'SE', 200, False, 'mappability_dedup_cds_RNA_SE_L200bp', 'transcript_uniqueness_factors_dedup_cds_RNA_SE_L200bp.tsv'),
    'v2_o': ('mane_rna',   'PE', 200, False, 'mappability_dedup_cds_RNA_PE_L200bp', 'transcript_uniqueness_factors_dedup_cds_RNA_PE_L200bp.tsv'),
}
V2_ORDER = ['v2_a', 'v2_b', 'v2_c', 'v2_d', 'v2_e', 'v2_f', 'v2_g', 'v2_h',
            'v2_i', 'v2_j', 'v2_k', 'v2_l', 'v2_m', 'v2_n', 'v2_o']


def section_header_v2(sid, src, mode, L, align_intron_max1, outdir_name):
    src_l = SRC_LABELS[src]
    return f'''#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# In silico mappability correction — {sid} — NF REFERENCE (MANE-restricted dedup CDS)
#
# Reads     : {src_l}
# Reference : dedup CDS, MANE_Select-restricted (nf) — NOT the genome-wide
#             Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna built for v7_nf/v8_b_nf
# Mode      : {mode}
# Read len  : {L} bp
# alignIntronMax=1 : {align_intron_max1}
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
RESULTS="$SCRIPT_DIR/../results"

REF="$SCRIPT_DIR/../../Ref"
GENOME_FA_NF="$REF/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
GTF_NF="$REF/Homo_sapiens.GRCh38.116.gtf"
MANE_RNA_NF="$REF/MANE_RNA_nf.fna"
DEDUP_CDS_MANE_NF_RAW="$REF/Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane_raw.fna"
DEDUP_CDS_MANE_NF="$REF/Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna"
CDS_MANE_IDX_NF="$REF/STAR_index_dedup_cds_mane_nf"

OUTDIR="$RESULTS/{outdir_name}_nf"
THREADS=64
READ_LEN={L}
FLANK=50

mkdir -p "$OUTDIR/simreads" "$OUTDIR/bam" "$OUTDIR/logs"

module purge
module load bio/STAR/2.7.11b-GCC-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0
module load lang/Python/3.11.5-GCCcore-13.2.0

ulimit -n 65536

echo "════════════════════════════════════════════════════════"
echo "Mappability — {sid} (nf reference, MANE-restricted dedup CDS, {src_l}, {mode}, {L} bp)"
echo "READ_LEN=$READ_LEN  THREADS=$THREADS"
echo "════════════════════════════════════════════════════════"
'''


def section_build_dedup_cds_mane_nf():
    """Step 0d(v2) — splice CDS (not full exon/UTR) from MANE_Select-tagged
    transcripts, then seqkit-dedup. Same resumable local-seqkit fallback as
    the genome-wide Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds.fna build. Shared across all v2_nf scripts;
    skipped if Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna already exists."""
    return r'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 0d: Build Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna — MANE_Select-restricted CDS,
# deduplicated with seqkit (shared; skipped if already built)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0d: Checking Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna..."
if [ -f "$DEDUP_CDS_MANE_NF" ]; then
    echo "Step 0d: Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna already exists, skipping."
    SKIP_STEP_0D=1
else
    SKIP_STEP_0D=0
fi

if [ "$SKIP_STEP_0D" = "0" ] && [ -f "$DEDUP_CDS_MANE_NF_RAW" ]; then
    echo "Step 0d: Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane_raw.fna already extracted, skipping extraction."
elif [ "$SKIP_STEP_0D" = "0" ]; then
    echo "Step 0d: Extracting MANE_Select CDS from nf GTF..."
python3 - << PYEOF
import pysam
from pathlib import Path
from collections import defaultdict

gtf_path  = Path("$GTF_NF")
genome_fa = Path("$GENOME_FA_NF")
raw_fa    = Path("$DEDUP_CDS_MANE_NF_RAW")

def parse_gtf_attrs(col9):
    attrs = defaultdict(list)
    for field in col9.strip().split(';'):
        field = field.strip()
        if not field:
            continue
        key, _, val = field.partition(' ')
        attrs[key].append(val.strip('"'))
    return attrs

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def revcomp(s): return s.translate(RC)[::-1]

mane_select_tx = set()
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'transcript':
            continue
        a = parse_gtf_attrs(f[8])
        if 'MANE_Select' not in a.get('tag', []):
            continue
        tid = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None:
            continue
        mane_select_tx.add(f"{tid}.{tver}" if tver else tid)
print(f"  MANE_Select transcripts: {len(mane_select_tx):,}", flush=True)

tx_cds = defaultdict(list)
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        a = parse_gtf_attrs(f[8])
        tid = a.get('transcript_id', [None])[0]
        tver = a.get('transcript_version', [None])[0]
        if tid is None:
            continue
        versioned_tid = f"{tid}.{tver}" if tver else tid
        if versioned_tid not in mane_select_tx:
            continue
        chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
        tx_cds[versioned_tid].append((chrom, start, end, strand))

for tid in tx_cds:
    tx_cds[tid].sort(key=lambda e: e[1])
print(f"  MANE_Select transcripts with CDS: {len(tx_cds):,}", flush=True)

fa = pysam.FastaFile(str(genome_fa))
n_written = n_skip = 0
with open(raw_fa, 'w') as out:
    for tid, segs in tx_cds.items():
        strand = segs[0][3]
        try:
            seq = ''.join(fa.fetch(c, s, e) for (c, s, e, _st) in segs).upper()
        except (KeyError, ValueError):
            n_skip += 1
            continue
        if not seq:
            n_skip += 1
            continue
        if strand == '-':
            seq = revcomp(seq)
        out.write(f">{tid}\n{seq}\n")
        n_written += 1
fa.close()

print(f"  CDS sequences written (pre-dedup): {n_written:,}  (skipped: {n_skip})")
print(f"  Saved raw: {raw_fa}", flush=True)
PYEOF
fi

if [ "$SKIP_STEP_0D" = "0" ]; then
    if command -v seqkit &> /dev/null; then
        echo "  Deduplicating exact-sequence duplicates with seqkit rmdup -s..."
        seqkit rmdup -s -i "$DEDUP_CDS_MANE_NF_RAW" -o "$DEDUP_CDS_MANE_NF" --threads $THREADS
        echo "Step 0d complete."
    else
        echo ""
        echo "  seqkit not found on PATH on this cluster."
        echo "  $DEDUP_CDS_MANE_NF_RAW is ready (MANE_Select CDS, not yet deduplicated)."
        echo "  To finish Step 0d:"
        echo "    1. scp $DEDUP_CDS_MANE_NF_RAW to a machine with seqkit installed"
        echo "    2. run: seqkit rmdup -s -i Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane_raw.fna -o Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna"
        echo "    3. scp Homo_sapiens.GRCh38.dna_sm.primary_ass.dedup_cds_mane.fna back to: $DEDUP_CDS_MANE_NF"
        echo "    4. re-run this script — it will detect the file and skip straight past Step 0d"
        exit 1
    fi
fi
'''


def section_build_cds_mane_index_nf():
    """Step 0e(v2) — build STAR_index_dedup_cds_mane_nf (shared; skipped if exists)."""
    return r'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 0e: Build STAR index (dedup CDS, MANE-restricted, nf) — shared; skipped if exists
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 0e: Checking STAR dedup-CDS-MANE-nf index..."
mkdir -p "$CDS_MANE_IDX_NF"
if [ ! -f "$CDS_MANE_IDX_NF/SA" ]; then
    STAR \
        --runMode             genomeGenerate \
        --genomeDir           "$CDS_MANE_IDX_NF" \
        --genomeFastaFiles    "$DEDUP_CDS_MANE_NF" \
        --genomeSAindexNbases 12 \
        --genomeChrBinNbits   11 \
        --runThreadN          $THREADS \
        2>&1 | tee "$OUTDIR/logs/star_index_cds_mane_nf.log"
    echo "Step 0e complete."
else
    echo "Step 0e: index already exists, skipping."
fi
'''


def section_star_align_v2(mode, align_intron_max1):
    """Step 2(v2) — STAR alignment against the MANE-restricted dedup-CDS-nf index."""
    extra = '    --alignIntronMax        1 \\\n' if align_intron_max1 else ''
    if mode == 'PE':
        reads_in = (
            '    --readFilesIn           "$OUTDIR/simreads/sim_R1.fastq" \\\n'
            '                            "$OUTDIR/simreads/sim_R2.fastq" \\'
        )
    else:
        reads_in = '    --readFilesIn           "$OUTDIR/simreads/sim_R1.fastq" \\'

    return f'''
# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: STAR alignment — {mode}, MANE-restricted dedup CDS (nf)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Step 2: STAR alignment ({mode}, dedup-CDS-MANE-nf, $THREADS threads)..."

STAR \\
    --runThreadN            $THREADS \\
    --genomeDir             "$CDS_MANE_IDX_NF" \\
{reads_in}
    --outSAMtype            BAM SortedByCoordinate \\
    --outSAMattributes      NH HI AS NM \\
    --outSAMmultNmax        1 \\
    --outFilterMultimapNmax 20 \\
{extra}    --outBAMsortingThreadN  $THREADS \\
    --outBAMsortingBinsN    20 \\
    --limitBAMsortRAM       160000000000 \\
    --outFileNamePrefix     "$OUTDIR/bam/sim_" \\
    2>&1 | tee "$OUTDIR/logs/star_sim.log"

samtools index -@ $THREADS "$OUTDIR/bam/sim_Aligned.sortedByCoord.out.bam"
echo "Step 2 complete."
'''


def make_script_v2(sid):
    src, mode, L, align_intron_max1, outdir_name, old_tsv = OLD_V2[sid]
    tsv_n = old_tsv.replace('.tsv', '_nf.tsv')
    parts = [section_header_v2(sid, src, mode, L, align_intron_max1, outdir_name)]
    if src == 'exon_flank':
        parts.append(section_chr_map())
    parts.append(section_build_dedup_cds_mane_nf())
    parts.append(section_build_cds_mane_index_nf())
    if src == 'mane_rna':
        parts.append(section_build_mane_rna())
    if src == 'exon_flank':
        parts.append(section_simulate_exon_flank(mode))
    else:
        parts.append(section_simulate_mane_rna(mode))
    parts.append(section_star_align_v2(mode, align_intron_max1))
    parts.append(section_compute_uniqueness(mode, tsv_n))
    return ''.join(parts)


if __name__ == '__main__':
    created = []
    for sid, src, mode, L in NEW_SCRIPTS_NF:
        content  = make_script(sid, src, mode, L)
        out_path = HERE / f'mappability_correction_{sid}.sh'
        out_path.write_text(content)
        out_path.chmod(out_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        created.append(sid)
        print(f"  {sid:10s}  {SRC_LABELS[src]:<20}  {mode}  {L:>3} bp")

    v8_content = make_script_v8a_nf()
    v8_path = HERE / 'mappability_correction_v8_a_nf.sh'
    v8_path.write_text(v8_content)
    v8_path.chmod(v8_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    created.append('v8_a_nf')
    print(f"  {'v8_a_nf':10s}  {SRC_LABELS['mane_rna_proper']:<20}  PE  150 bp (proper)")

    for sid, mode, L, target in NEW_SCRIPTS_V7_NF:
        content  = make_script_v7(sid, mode, L, target)
        out_path = HERE / f'mappability_correction_{sid}.sh'
        out_path.write_text(content)
        out_path.chmod(out_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        created.append(sid)
        print(f"  {sid:10s}  dedup CDS reads      {mode}  {L:>3} bp  ({target})")

    for which, sid in (('b', 'v8_b_nf'), ('c', 'v8_c_nf')):
        content  = make_script_v8bc_nf(which)
        out_path = HERE / f'mappability_correction_{sid}.sh'
        out_path.write_text(content)
        out_path.chmod(out_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        created.append(sid)
        print(f"  {sid:10s}  MANE RNA (proper PE) reuse from v8_a_nf")

    for sid in V2_ORDER:
        nf_sid = f'{sid}_nf'
        content  = make_script_v2(sid)
        out_path = HERE / f'mappability_correction_{nf_sid}.sh'
        out_path.write_text(content)
        out_path.chmod(out_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        created.append(nf_sid)
        src, mode, L, aim1, _, _ = OLD_V2[sid]
        print(f"  {nf_sid:10s}  {SRC_LABELS[src]:<20}  {mode}  {L:>3} bp  (dedup CDS, MANE-restricted{' , alignIntronMax1' if aim1 else ''})")

    print(f"\nCreated {len(created)} scripts.")
