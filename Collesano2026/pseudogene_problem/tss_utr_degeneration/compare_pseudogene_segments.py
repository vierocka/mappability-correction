#!/usr/bin/env python3
"""Does a chr6 MHC pseudogene keep its parent's CDS but lose the promoter/TSS and UTRs?

For each pseudogene, take the best-matching functional parent among the candidate genes. Split the parent's
MANE Select transcript (GENCODE v49) into promoter (1 kb upstream of the TSS), 5'UTR, CDS and 3'UTR.
Align each segment (unmasked, upper-cased) with blastn to the pseudogene locus ±5 kb on the primary assembly.
conserved = mean over segment bases of the identity of the best HSP covering the base (0 if uncovered).

usage: compare_pseudogene_segments.py GENCODE_GTF.gz GENOME.fa(.fai) OUTDIR
figure: plot_pseudogene_segments.py OUTDIR
"""
import gzip, re, subprocess, sys, tempfile
from pathlib import Path

GTF, GENOME, OUT = sys.argv[1], sys.argv[2], Path(sys.argv[3])
PAIRS = {  # pseudogene: candidate parents (best CDS match is used)
    "HLA-H": ["HLA-A", "HLA-B", "HLA-C"], "HLA-J": ["HLA-A", "HLA-B", "HLA-C"],
    "HLA-K": ["HLA-A", "HLA-B", "HLA-C"], "HLA-L": ["HLA-A", "HLA-B", "HLA-C"],
    "HLA-W": ["HLA-A", "HLA-B", "HLA-C"], "HLA-DRB6": ["HLA-DRB1"],
    "HLA-DPB2": ["HLA-DPB1"], "HLA-DPA2": ["HLA-DPA1"], "HLA-DQB3": ["HLA-DQB1"],
    "MICE": ["MICA", "MICB"],
}
FLANK, PROM = 5000, 1000
need = set(PAIRS) | {p for v in PAIRS.values() for p in v}
genes, tx = {}, {}
attr = lambda s, k: (re.search(k + r' "([^"]+)"', s) or [None, None])[1]
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if c[0] != "chr6": continue
        g = attr(c[8], "gene_name")
        if g not in need: continue
        s, e, st = int(c[3]), int(c[4]), c[6]
        if c[2] == "gene":
            genes[g] = (s, e, st, attr(c[8], "gene_type"))
        elif 'tag "MANE_Select"' in c[8] and c[2] in ("exon", "CDS", "UTR", "start_codon", "stop_codon"):
            tx.setdefault(g, {"st": st, "exon": [], "CDS": [], "stop": []})
            if c[2] in ("exon", "CDS"): tx[g][c[2]].append((s, e))
            if c[2] == "stop_codon": tx[g]["stop"].append((s, e))

def fetch(region):
    r = subprocess.run(["samtools", "faidx", GENOME, region], capture_output=True, text=True, check=True)
    return "".join(r.stdout.split("\n")[1:])

rc = lambda s: s[::-1].translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))

def segments(g):
    t = tx[g]; st = t["st"]
    ex = sorted(t["exon"]); cds = sorted(t["CDS"] + t["stop"])
    cs, ce = min(a for a, b in cds), max(b for a, b in cds)
    def spliced(lo, hi):
        parts = [fetch(f"6:{max(a, lo)}-{min(b, hi)}") for a, b in ex if min(b, hi) >= max(a, lo)]
        seq = "".join(parts); return rc(seq) if st == "-" else seq
    first, last = ex[0][0], ex[-1][1]
    if st == "+":
        prom = fetch(f"6:{first - PROM}-{first - 1}")
        u5, u3 = spliced(first, cs - 1), spliced(ce + 1, last)
    else:
        prom = rc(fetch(f"6:{last + 1}-{last + PROM}"))
        u5, u3 = spliced(ce + 1, last), spliced(first, cs - 1)
    return {"promoter": prom, "5'UTR": u5, "CDS": spliced(cs, ce), "3'UTR": u3}

def conserved(q, target_fa, tmp):
    if len(q) < 50: return float("nan")      # too short to align reliably
    qf = Path(tmp) / "q.fa"; qf.write_text(">q\n" + q + "\n")
    r = subprocess.run(["blastn", "-task", "blastn", "-query", str(qf), "-subject", str(target_fa),
                        "-evalue", "1e-5", "-outfmt", "6 qstart qend pident bitscore"],
                       capture_output=True, text=True, check=True)
    hsps = sorted((tuple(map(float, l.split())) for l in r.stdout.splitlines()), key=lambda h: -h[3])
    score = [0.0] * len(q)           # each query base scored once, best HSP first
    for qs, qe, pid, _ in hsps:
        for i in range(int(min(qs, qe)) - 1, int(max(qs, qe))):
            if score[i] == 0.0: score[i] = pid / 100
    return sum(score) / len(q)

OUT.mkdir(parents=True, exist_ok=True)
rows = []
with tempfile.TemporaryDirectory() as tmp:
    for pg, cands in PAIRS.items():
        s, e, st, gt = genes[pg]
        tfa = Path(tmp) / "t.fa"; tfa.write_text(">t\n" + fetch(f"6:{s - FLANK}-{e + FLANK}") + "\n")
        best = max(cands, key=lambda p: conserved(segments(p)["CDS"].upper(), tfa, tmp))
        seg = segments(best)
        vals = {k: conserved(v.upper(), tfa, tmp) for k, v in seg.items()}
        rows.append([pg, gt, best] + [round(vals[k], 3) for k in ("promoter", "5'UTR", "CDS", "3'UTR")]
                    + [len(seg[k]) for k in ("5'UTR", "CDS", "3'UTR")])
        print(rows[-1], flush=True)

hdr = ["pseudogene", "gene_type", "parent", "promoter_1kb", "UTR5", "CDS", "UTR3", "len_UTR5", "len_CDS", "len_UTR3"]
(OUT / "pseudogene_segment_conservation.tsv").write_text(
    "\t".join(hdr) + "\n" + "\n".join("\t".join(map(str, r)) for r in rows) + "\n")

