#!/usr/bin/env python3
"""Ties between the top two kallisto isoforms per gene, and kallisto-max vs featureCounts agreement (HCC1937).
usage: isoform_ties.py abundance.tsv CDS_GRCh38p14_annotNames.csv fc_vs_kallisto_max_merged.tsv
Result 2026-11-22: 9,803 genes (>=2 isoforms with counts, >=10 counts); 30.7 % top-2 within 1 %;
Spearman 0.811; 735 genes fc=0 & kallisto>=10; 52 the reverse; median kallisto/fc 0.26."""
import sys, numpy as np, pandas as pd
from scipy.stats import spearmanr
ab, annot, merged = sys.argv[1:4]
a = pd.read_csv(annot, sep="\t").drop_duplicates("CDS").set_index("CDS").gene
k = pd.read_csv(ab, sep="\t")
k["iso"] = k.target_id.str.extract(r"_cds_(.+)_\d+$")[0]
k["gene"] = k.iso.map(a)
print("targets", len(k), "| targets sharing an accession (patch/alt copies):", k.duplicated("iso", keep=False).sum())
d = k.groupby(["gene", "iso"]).est_counts.max().reset_index()          # one value per accession
e = d[d.est_counts > 0].groupby("gene").est_counts.apply(lambda s: sorted(s, reverse=True))
e = e[e.apply(lambda s: sum(s) >= 10 and len(s) >= 2)]
tie = e.apply(lambda s: abs(s[0] - s[1]) / s[0] < 0.01)
print(f"genes: {len(e)}  top-2 tied within 1 %: {100 * tie.mean():.1f} %")
m = pd.read_csv(merged, sep="\t"); x = m[(m.fc_counts > 0) | (m.kallisto_max_est_counts > 0)]
b = (x.fc_counts >= 10) & (x.kallisto_max_est_counts >= 10)
print(f"Spearman {spearmanr(x.fc_counts, x.kallisto_max_est_counts)[0]:.3f}",
      "| fc=0 & kal>=10:", ((x.fc_counts == 0) & (x.kallisto_max_est_counts >= 10)).sum(),
      "| kal=0 & fc>=10:", ((x.kallisto_max_est_counts == 0) & (x.fc_counts >= 10)).sum(),
      f"| median kal/fc {np.median(x.kallisto_max_est_counts[b] / x.fc_counts[b]):.2f}")
