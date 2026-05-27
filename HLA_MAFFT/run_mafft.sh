#!/bin/bash
# Align each HLA FASTA with MAFFT (--auto mode, 4 threads)
set -euo pipefail

OUT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for fna in \
    "$OUT/HLA_MANE.fna" \
    "$OUT/HLA_dedup_cds.fna" \
    "$OUT/HLA_primaryChr.fna"; do

    base="${fna%.fna}"
    name="$(basename "$fna")"
    echo "━━━━ MAFFT: $name"
    mafft --auto --thread 4 "$fna" \
        > "${base}.mafft.fna" \
        2> "${base}.mafft.log"
    n=$(grep -c '^>' "${base}.mafft.fna")
    echo "  Sequences aligned: $n"
    echo "  → ${base}.mafft.fna"
    echo ""
done

echo "All alignments complete."
