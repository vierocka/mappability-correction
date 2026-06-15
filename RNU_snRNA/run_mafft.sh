#!/usr/bin/env bash
# Per-family MAFFT alignment with high-accuracy settings.
# --globalpair --maxiterate 1000 = L-INS-i (global alignment, iterative refinement).
# Appropriate for short (63–191 bp), highly similar sequences within each family.
# Cross-family alignment is deliberately avoided — it only produces uninformative gaps.
set -euo pipefail

DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"

FAMILIES=(RNU1 RNU2 RNU4 RNU5 RNU6)

for FAM in "${FAMILIES[@]}"; do
    IN="${FAM}_family.fna"
    OUT="${FAM}_family.mafft.fna"
    LOG="${FAM}_family.mafft.log"
    echo "=== MAFFT: ${FAM} family ==="
    mafft --globalpair --maxiterate 1000 --quiet "$IN" > "$OUT" 2> "$LOG"
    echo "  Done → $OUT"
done

echo ""
echo "All per-family alignments complete."
