#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=4:00:00
#SBATCH --output=mappability_v7_fix_uf_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL

# ─────────────────────────────────────────────────────────────────────────────
# Fix per-CDS uniqueness factor computation for all v7 variants.
#
# The original scripts grouped all reads under cds_id='lcl' because
#   line[1:].split('|')[0]  →  'lcl'  (from 'lcl|NC_..._NP_...|pos')
# The fix uses
#   '|'.join(line[1:].split('|')[:-1])  →  'lcl|NC_..._NP_...'
#
# Re-runs Step 3 only — BAM files and simreads already exist.
# No STAR re-alignment needed.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

BASE="./mappability-correction"
RESULTS="$BASE/mappability/results"

module purge
module load lang/Python/3.11.5-GCCcore-13.2.0

compute_uf_fixed () {
    local OUTDIR="$1"
    local TSV_NAME="$2"
    local LABEL="$3"

    local BAM="$OUTDIR/bam/sim_Aligned.sortedByCoord.out.bam"
    local R1="$OUTDIR/simreads/sim_R1.fastq"

    if [ ! -f "$BAM" ]; then
        echo "  SKIP $LABEL — BAM not found: $BAM"
        return
    fi
    if [ ! -f "$R1" ]; then
        echo "  SKIP $LABEL — reads not found: $R1"
        return
    fi

    echo "  Computing: $LABEL"

python3 - << PYEOF
import pysam, csv
from collections import defaultdict
from pathlib import Path

outdir   = Path("$OUTDIR")
bam_path = outdir / "bam/sim_Aligned.sortedByCoord.out.bam"
r1_path  = outdir / "simreads/sim_R1.fastq"

# Count simulated reads per CDS sequence.
# Read name format: lcl|NC_..._NP_...|pos  — strip trailing |pos
sim_counts = defaultdict(int)
with open(r1_path) as fh:
    for line in fh:
        if line.startswith('@'):
            parts = line[1:].rstrip().split('|')
            cid = '|'.join(parts[:-1])
            sim_counts[cid] += 1

print(f"  CDS sequences: {len(sim_counts):,}   R1 reads: {sum(sim_counts.values()):,}")

# Count uniquely mapping reads per CDS sequence from BAM.
unique_back = defaultdict(int)
total_back  = defaultdict(int)
bam = pysam.AlignmentFile(str(bam_path), 'rb')
for read in bam.fetch():
    if read.is_supplementary or read.is_secondary or read.is_unmapped:
        continue
    if read.is_read2:
        continue
    parts = read.query_name.split('|')
    cid   = '|'.join(parts[:-1])
    total_back[cid] += 1
    if read.mapping_quality == 255:
        unique_back[cid] += 1
bam.close()

records = []
for cid, n_sim in sim_counts.items():
    n_u = unique_back.get(cid, 0)
    n_t = total_back.get(cid, 0)
    records.append({
        'cds_id':            cid,
        'n_positions':       n_sim,
        'n_unique_back':     n_u,
        'n_multi_back':      n_t - n_u,
        'n_unmapped':        n_sim - n_t,
        'uniqueness_factor': round(n_u / n_sim, 6) if n_sim > 0 else 0.0,
    })
records.sort(key=lambda r: r['uniqueness_factor'])

out = outdir / "$TSV_NAME"
fields = ['cds_id','n_positions','n_unique_back','n_multi_back','n_unmapped','uniqueness_factor']
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader()
    w.writerows(records)

fs = [r['uniqueness_factor'] for r in records]
n  = len(fs)
print(f"  n_cds={n:,}  mean_uf={sum(fs)/n:.4f}  "
      f"uf=1.0: {sum(1 for f in fs if f==1.0):,}  "
      f"uf=0.0: {sum(1 for f in fs if f==0.0):,}")
print(f"  Saved: {out}", flush=True)
PYEOF
}

echo "════════════════════════════════════════════════════════"
echo "v7 — fix per-CDS uniqueness factor computation"
echo "Started : $(date)"
echo "Host    : $(hostname)"
echo "════════════════════════════════════════════════════════"

echo ""
echo "── Full genome targets ──────────────────────────────────"
compute_uf_fixed "$RESULTS/mappability_dedup_cds_reads_SE"        cds_uniqueness_factors_genome_SE.tsv        "v7_a SE 100bp"
compute_uf_fixed "$RESULTS/mappability_dedup_cds_reads_PE"        cds_uniqueness_factors_genome_PE.tsv        "v7_b PE 100bp"
compute_uf_fixed "$RESULTS/mappability_dedup_cds_reads_SE_L75bp"  cds_uniqueness_factors_genome_SE_L75bp.tsv  "v7_c SE 75bp"
compute_uf_fixed "$RESULTS/mappability_dedup_cds_reads_PE_L75bp"  cds_uniqueness_factors_genome_PE_L75bp.tsv  "v7_d PE 75bp"
compute_uf_fixed "$RESULTS/mappability_dedup_cds_reads_SE_L150bp" cds_uniqueness_factors_genome_SE_L150bp.tsv "v7_e SE 150bp"
compute_uf_fixed "$RESULTS/mappability_dedup_cds_reads_PE_L150bp" cds_uniqueness_factors_genome_PE_L150bp.tsv "v7_f PE 150bp"
compute_uf_fixed "$RESULTS/mappability_dedup_cds_reads_SE_L200bp" cds_uniqueness_factors_genome_SE_L200bp.tsv "v7_i SE 200bp"
compute_uf_fixed "$RESULTS/mappability_dedup_cds_reads_PE_L200bp" cds_uniqueness_factors_genome_PE_L200bp.tsv "v7_j PE 200bp"

echo ""
echo "── CDS self-map targets ─────────────────────────────────"
compute_uf_fixed "$RESULTS/mappability_dedup_cds_reads_selfmap_SE" cds_uniqueness_factors_selfmap_SE.tsv "v7_g SE 100bp selfmap"
compute_uf_fixed "$RESULTS/mappability_dedup_cds_reads_selfmap_PE" cds_uniqueness_factors_selfmap_PE.tsv "v7_h PE 100bp selfmap"

echo ""
echo "════════════════════════════════════════════════════════"
echo "Done: $(date)"
echo ""
echo "Fixed TSVs:"
find "$RESULTS" -name "cds_uniqueness_factors_*.tsv" | sort
echo "════════════════════════════════════════════════════════"
