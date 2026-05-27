#!/bin/bash
# Align each HLA FASTA with MAFFT (--auto mode)
set -euo pipefail

OUT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THREADS="${1:-4}"   # pass thread count as first argument, default 4

for fna in \
    "$OUT/HLA_MANE.fna" \
    "$OUT/HLA_dedup_cds.fna" \
    "$OUT/HLA_primaryChr.fna"; do

    base="${fna%.fna}"
    name="$(basename "$fna")"
    echo "━━━━ MAFFT: $name"
    mafft --auto --thread "$THREADS" "$fna" \
        > "${base}.mafft.fna" \
        2> "${base}.mafft.log"
    n=$(grep -c '^>' "${base}.mafft.fna")
    echo "  Sequences aligned: $n"
    echo "  → ${base}.mafft.fna"
    echo ""
done

echo "All alignments complete."
