#!/usr/bin/env python3
"""Group transcripts whose terminal window is identical.

Transcripts sharing an identical last (or first) W bp cannot be told apart by any read falling
in that window. Grouping is by exact string match on the window, keyed by md5 - not an identity
threshold. A 3' window is the default because 3'-biased protocols (10x and similar) sequence it,
and because 3' UTRs are where paralogs usually diverge; --end 5 gives the 5' background.

Needs only a transcript FASTA, so it works for any organism. The transcript -> gene mapping is
taken from the FASTA header if present (GENCODE pipe format, or `gene=` / `gene_id=` fields),
otherwise from --tx2gene.

usage: terminal_clusters.py --fasta F --out F [--window 250] [--end 3] [--tx2gene F]
                            [--gene-table F] [--biotypes F]
"""
import argparse
import gzip
import re
from collections import defaultdict
from hashlib import md5

OPEN = lambda p: gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def read_fasta(path):
    name, chunks = None, []
    with OPEN(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def parse_header(h):
    """Return (transcript_id, gene_id, gene_name, biotype) from common header styles."""
    f = h.split("|")
    if len(f) >= 6 and f[0].startswith(("ENST", "ENSMUST", "ENS")):
        return f[0], f[1], f[5], (f[7] if len(f) > 7 else "")
    tid = h.split()[0]
    gid = (re.search(r'gene_id[= ]([^\s;]+)', h) or re.search(r'gene[= ]([^\s;]+)', h))
    bt = re.search(r'(?:transcript_biotype|gene_biotype|transcript_type)[= ]([^\s;]+)', h)
    gn = re.search(r'gene_symbol[= ]([^\s;]+)', h)
    return tid, (gid.group(1) if gid else ""), (gn.group(1) if gn else ""), (bt.group(1) if bt else "")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--gene-table", help="also write per-gene cluster summary here")
    p.add_argument("--window", type=int, default=250)
    p.add_argument("--end", choices=["3", "5"], default="3")
    p.add_argument("--tx2gene", help="TSV transcript_id<TAB>gene_id[<TAB>biotype]")
    a = p.parse_args()

    ext = {}
    if a.tx2gene:
        with open(a.tx2gene) as fh:
            for line in fh:
                if line.startswith(("transcript", "#")):
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) >= 2:
                    ext[c[0]] = (c[1], c[2] if len(c) > 2 else "")

    win_to_tx = defaultdict(list)
    info = {}
    n_short = 0
    for h, seq in read_fasta(a.fasta):
        tid, gid, gname, bt = parse_header(h)
        if tid in ext:
            gid = ext[tid][0] or gid
            bt = ext[tid][1] or bt
        # A transcript shorter than the window is keyed on its whole sequence rather than
        # dropped, so identical short duplicates (yeast CUP1, the ASP3 array) still cluster.
        # Keys of different length never collide, so this adds no false pairs. The residual
        # blind spot: a short transcript identical to the END of a longer one is not detected.
        if len(seq) < a.window:
            n_short += 1
        w = seq[-a.window:] if a.end == "3" else seq[:a.window]
        win_to_tx[w].append(tid)
        info[tid] = (gid or tid, gname or gid or tid, bt)

    if not info:
        raise SystemExit(f"no usable sequence in {a.fasta} (every record shorter than the "
                         f"{a.window} bp window, or the file is empty)")
    clusters = sorted(((w, t) for w, t in win_to_tx.items() if len(t) > 1), key=lambda x: -len(x[1]))

    is_pseudo = lambda t: "pseudogene" in (t or "")
    rows, gene_rows = [], defaultdict(lambda: {"cluster_ids": set(), "n_tx": 0,
                                               "n_pseudo_partners": 0, "partner_genes": set()})
    for i, (w, tids) in enumerate(clusters, 1):
        cid = f"C{i}_{md5(w.encode()).hexdigest()[:10]}"
        genes = {info[t][0] for t in tids}
        bts = {info[t][2] for t in tids}
        mixed = int(any(is_pseudo(b) for b in bts) and any(not is_pseudo(b) and b for b in bts))
        for t in tids:
            gid, gname, bt = info[t]
            rows.append([cid, len(tids), len(genes), t, gid, gname, bt, mixed])
            g = gene_rows[gid]
            g["cluster_ids"].add(cid)
            g["n_tx"] += 1
            g["n_pseudo_partners"] += sum(is_pseudo(info[o][2]) and info[o][0] != gid for o in tids)
            g["partner_genes"] |= {info[o][1] for o in tids if info[o][0] != gid}

    with open(a.out, "w") as fh:
        fh.write("cluster_id\tcluster_size\tn_genes\ttranscript_id\tgene_id\tgene_name\tbiotype\tmixed_coding_pseudogene\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    if a.gene_table:
        with open(a.gene_table, "w") as fh:
            fh.write("gene_id\tn_clusters\tn_transcripts_in_clusters\tn_pseudogene_partners\tn_partner_genes\tpartner_genes\n")
            for gid, g in sorted(gene_rows.items()):
                fh.write(f"{gid}\t{len(g['cluster_ids'])}\t{g['n_tx']}\t{g['n_pseudo_partners']}\t"
                         f"{len(g['partner_genes'])}\t{','.join(sorted(g['partner_genes'])[:20])}\n")

    pct = 100 * n_short / max(len(info), 1)
    print(f"{len(clusters)} clusters of identical {a.window} bp {a.end}' ends "
          f"covering {len(rows)} transcripts and {len(gene_rows)} genes")
    print(f"{n_short} transcripts ({pct:.1f} %) are shorter than the {a.window} bp window and "
          f"were keyed on their whole sequence")
    if pct > 5:
        print(f"NOTE: {pct:.1f} % of transcripts fall below the window. A short transcript that "
              f"matches the END of a longer one is not detected, so causes may be "
              f"under-determined for them. Consider a smaller --window for this annotation.")


if __name__ == "__main__":
    main()
