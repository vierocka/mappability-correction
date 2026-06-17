#!/bin/bash
# Download cell-type marker tables for 3' UF impact analysis.
#
# Files are already present in this directory (downloaded June 2026).
# This script re-downloads them if missing or if you want to update.
#
# Sources:
#   CellMarker 2.0  — Chen et al., NAR 2023, PMID 36300619
#     http://bio-bigdata.hrbmu.edu.cn/CellMarker/
#   PanglaoDB       — Franzén et al., Database 2019, PMID 30951143
#     https://panglaodb.se/

set -euo pipefail
cd "$(dirname "$0")"

echo "Downloading CellMarker 2.0 (human cell markers XLSX)..."
curl -fL -o CellMarker2_Human_cell_markers.xlsx \
    "http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download_files/file/Cell_marker_Human.xlsx"

echo "Downloading PanglaoDB markers (2020-03-27 release)..."
curl -fL -o PanglaoDB_markers_27Mar2020.tsv.gz \
    "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz" \
  && gunzip -f PanglaoDB_markers_27Mar2020.tsv.gz

echo ""
echo "Checking files:"
python3 - << 'PYEOF'
import openpyxl, gzip, sys

xlsx = "CellMarker2_Human_cell_markers.xlsx"
wb = openpyxl.load_workbook(xlsx)
ws = wb.active
n = sum(1 for row in ws.iter_rows(min_row=2, values_only=True) if any(v is not None for v in row))
print(f"CellMarker2 XLSX  : {n:,} data rows, columns: {[c.value for c in ws[1]][:5]}...")

with gzip.open("PanglaoDB_markers_27Mar2020.tsv.gz", "rt") as fh:
    header = fh.readline().rstrip().split("\t")
    n_rows = sum(1 for _ in fh)
print(f"PanglaoDB TSV.GZ  : {n_rows:,} data rows, columns: {header[:5]}...")
PYEOF

echo "Done."
