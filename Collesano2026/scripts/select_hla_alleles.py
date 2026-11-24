#!/usr/bin/env python3
"""
2-field HLA calls -> <sample>_alleles.fa (IPD-IMGT/HLA hla_nuc records) for build_personal_MANE_utr.py.
Representative per 2-field allele: lowest-numbered full allele (e.g. A*02:01:01:01) whose CDS
has the MANE CDS length (A 1098, B 1089, C 1101) and starts ATG / ends with a stop.
Same rule for every sample (old JY used C*07:02:01:03; this picks :01:01; identical protein).
Homozygous gene -> one record (MANE ID kept); heterozygous -> two (allele-suffixed IDs).
"""
import argparse
import re
import sys

CDS_LEN = {"A": 1098, "B": 1089, "C": 1101}
STOP = {"TAA", "TAG", "TGA"}


def read_nuc(path, genes):
    recs, h, s = {}, None, []
    for ln in open(path):
        if ln.startswith(">"):
            if h:
                recs[h] = "".join(s)
            h, s = ln[1:].strip(), []
        else:
            s.append(ln.strip())
    if h:
        recs[h] = "".join(s)
    out = {}
    for h, seq in recs.items():
        name = h.split()[1]
        if name.split("*")[0] in genes:
            out[name] = (h, seq.upper())
    return out


def key(name):
    f = re.sub(r"[A-Z]$", "", name.split("*")[1]).split(":")   # drop expression suffix for sorting
    return [int(x) if x.isdigit() else 999 for x in f], name


def pick(two_field, nuc):
    gene = two_field.split("*")[0]
    cands = [n for n in nuc if n == two_field or n.startswith(two_field + ":")]
    cands = [n for n in cands if not re.search(r"[NQLSC]$", n)]               # no null/questionable
    ok = [n for n in cands if len(nuc[n][1]) == CDS_LEN[gene]
          and nuc[n][1][:3] == "ATG" and nuc[n][1][-3:] in STOP]
    if not ok:
        raise SystemExit(f"{two_field}: no full-length CDS in hla_nuc ({len(cands)} candidates)")
    return sorted(ok, key=key)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True, help="final calls TSV from reconcile_hla.py (gene, allele1, allele2)")
    ap.add_argument("--hla-nuc", required=True, help="IPD-IMGT/HLA hla_nuc.fasta")
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", required=True)
    a = ap.parse_args()
    import pandas as pd
    calls = pd.read_csv(a.calls, sep="\t")
    nuc = read_nuc(a.hla_nuc, set(CDS_LEN))
    rows, fa = [], []
    for _, r in calls.iterrows():
        alle = [r.allele1] if r.allele1 == r.allele2 else [r.allele1, r.allele2]
        for t in alle:
            full = pick(t, nuc)
            h, seq = nuc[full]
            fa.append(f">{h}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)))
            rows.append({"gene": r.gene, "call": t, "imgt_allele": full, "imgt_id": h.split()[0],
                         "cds_bp": len(seq), "zygosity": "hom" if len(alle) == 1 else "het"})
    open(a.out, "w").write("\n".join(fa) + "\n")
    pd.DataFrame(rows).to_csv(a.report, sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
