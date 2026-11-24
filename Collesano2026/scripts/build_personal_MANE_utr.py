#!/usr/bin/env python3
"""
Build a personalized MANE FASTA by replacing the CDS of HLA-A/B/C with
sample-specific alleles, KEEPING the reference MANE 5'UTR and 3'UTR.

Supersedes {JY,HCC1937,HCT116}_build_personal_MANE.py (2026-10-05). Those
replaced the WHOLE MANE transcript with the bare IPD-IMGT/HLA CDS (hla_nuc,
ATG..stop), so every personalized HLA-A/B/C record lost its UTRs -- in ALL
three samples, JY included (homozygous entries kept the original ID/header
but still got the CDS-only sequence). Since cds_mane.saf counts every other
gene over its whole transcript (Start=1, End=transcript length -- the name
is historical), HLA reads falling in the UTRs were lost and HLA-A/B/C were
under-counted relative to every other gene (~29% of transcript length).

UTR source: the MANE reference transcript itself (same UTRs the base run
uses). MANE CDS boundaries are taken from the MANE GFF and verified here
(start codon ATG, stop codon, 5'UTR + CDS + 3'UTR == transcript length,
allele CDS length == MANE CDS length). Checked against IPD-IMGT/HLA 3.63.1
hla_gen.fasta (2026-10-05): allele-specific UTRs differ from MANE's by 0-3
nt (5'UTR) and 0-6 nt (3'UTR) for HLA-A/B and C*07, but by ~12 SNVs + a
1-nt indel for C*03:04 and C*05:01 (3'UTR, 420 nt). That is within STAR's
default mismatch tolerance, so MANE UTRs are used as-is for every allele.

Each HLA record gets a " cds:<start>-<end>" header tag (1-based, inclusive,
transcript coordinates) so downstream scripts can slice the allele CDS
back out (build_personal_MANE_protein.py, build_allele_mappability_simreads.py).
STAR and featureCounts only use the first header token, so IDs are unchanged.

When two alleles exist for a gene, both are included as separate entries
(e.g. ENST00000376809.10_A2301 / _A2402), exactly as before.

Usage:
    python3 build_personal_MANE_utr.py --sample HCC1937
    python3 build_personal_MANE_utr.py --sample HCT116
    python3 build_personal_MANE_utr.py --sample JY
"""

import argparse
import re
from pathlib import Path

HLA_TRANSCRIPTS = {
    "HLA-A": "ENST00000376809.10",
    "HLA-B": "ENST00000412585.7",
    "HLA-C": "ENST00000376228.10",
}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def parse_fasta(path):
    """Return list of (header, seq) tuples."""
    records = []
    header = seq_parts = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def write_fasta(records, path, line_width=60):
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + "\n")


def allele_gene(header):
    """IPD-IMGT/HLA header -> HLA-A / HLA-B / HLA-C."""
    m = re.search(r'\b([ABC])\*', header)
    return f"HLA-{m.group(1)}" if m else None


def allele_short_id(header):
    """'A*23:01:01:01' -> 'A2301'."""
    m = re.search(r'\b([ABC]\*\d+:\d+)', header)
    return re.sub(r'[^A-Za-z0-9]', '', m.group(1)) if m else None


def mane_gene_symbol(header):
    m = re.search(r'gene_symbol:(HLA-[ABC])\b', header)
    return m.group(1) if m else None


def utr_lengths_from_gff(gff_path, transcript_ids):
    """transcript_id -> (len 5'UTR, len CDS, len 3'UTR), summed over GFF parts."""
    lens = {t: [0, 0, 0] for t in transcript_ids}
    col = {"five_prime_UTR": 0, "CDS": 1, "three_prime_UTR": 2}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] not in col:
                continue
            m = re.search(r'transcript_id=([^;]+)', p[8])
            if not m:
                continue
            tid = m.group(1)
            hit = next((t for t in transcript_ids
                        if t == tid or t.split(".")[0] == tid.split(".")[0]), None)
            if hit:
                lens[hit][col[p[2]]] += int(p[4]) - int(p[3]) + 1
    return {t: tuple(v) for t, v in lens.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True, help="JY / HCC1937 / HCT116")
    ap.add_argument("--mane",    default="MANE.GRCh38.v1.5.ensembl_rna.fna")
    ap.add_argument("--gff",     default="MANE.GRCh38.v1.5.ensembl_genomic.gff")
    ap.add_argument("--saf",     default="cds_mane.saf")
    ap.add_argument("--alleles", default=None, help="default: <sample>_alleles.fa")
    ap.add_argument("--out",     default=None, help="default: <sample>_personal_MANE.fa")
    ap.add_argument("--out_saf", default=None, help="default: <sample>_cds_mane.saf")
    args = ap.parse_args()
    alleles_fa = args.alleles or f"{args.sample}_alleles.fa"
    out_fa     = args.out     or f"{args.sample}_personal_MANE.fa"
    out_saf    = args.out_saf or f"{args.sample}_cds_mane.saf"

    # --- load sample alleles (CDS only, IPD-IMGT/HLA hla_nuc) --------------
    allele_seqs = {}   # "HLA-A" -> [(allele_short_id, cds), ...]
    for header, seq in parse_fasta(alleles_fa):
        gene, aid = allele_gene(header), allele_short_id(header)
        if not (gene and aid):
            raise SystemExit(f"Cannot determine gene/allele for header: {header}")
        allele_seqs.setdefault(gene, []).append((aid, seq.upper()))
        print(f"  loaded allele {gene}: {header.split()[0]} {aid} ({len(seq)} bp)")

    # --- MANE UTR/CDS geometry from the GFF --------------------------------
    geometry = utr_lengths_from_gff(args.gff, set(HLA_TRANSCRIPTS.values()))

    mane_records = parse_fasta(args.mane)
    out_records = []
    saf_replacements = {}   # base transcript_id -> [(new_id, new_len), ...]

    for header, seq in mane_records:
        gene = mane_gene_symbol(header)
        if not (gene and gene in allele_seqs):
            out_records.append((header, seq))
            continue

        transcript_id = header.split()[0]
        rest = header[len(transcript_id):].rstrip()
        u5, c, u3 = geometry[transcript_id]
        # verify MANE geometry before trusting it
        if u5 + c + u3 != len(seq):
            raise SystemExit(f"{transcript_id}: 5'UTR {u5} + CDS {c} + 3'UTR {u3} "
                             f"!= transcript length {len(seq)}")
        ref_cds = seq[u5:u5 + c]
        if ref_cds[:3] != "ATG" or ref_cds[-3:] not in STOP_CODONS:
            raise SystemExit(f"{transcript_id}: MANE CDS slice does not start with ATG "
                             f"/ end with a stop ({ref_cds[:3]}..{ref_cds[-3:]})")
        utr5, utr3 = seq[:u5], seq[u5 + c:]

        alleles = allele_seqs[gene]
        replacement_list = []
        for aid, cds in alleles:
            if len(cds) != c or cds[:3] != "ATG" or cds[-3:] not in STOP_CODONS:
                raise SystemExit(f"{gene} {aid}: allele CDS ({len(cds)} bp, "
                                 f"{cds[:3]}..{cds[-3:]}) incompatible with MANE CDS ({c} bp)")
            new_seq = utr5 + cds + utr3
            cds_tag = f" cds:{u5 + 1}-{u5 + c}"
            if len(alleles) == 1:
                new_id = transcript_id                       # homozygous: ID unchanged
                new_header = f"{transcript_id}{rest}{cds_tag}"
            else:
                new_id = f"{transcript_id}_{aid}"
                new_header = f"{new_id}{rest} allele:{aid}{cds_tag}"
            out_records.append((new_header, new_seq))
            replacement_list.append((new_id, len(new_seq)))
            print(f"  {gene} {new_id}: 5'UTR {u5} + CDS {len(cds)} [{aid}] "
                  f"+ 3'UTR {u3} = {len(new_seq)} bp (MANE {len(seq)} bp)")
        saf_replacements[transcript_id] = replacement_list

    missing = set(allele_seqs) - {mane_gene_symbol(h) for h, _ in mane_records}
    if missing:
        raise SystemExit(f"No MANE transcript found for {missing}")

    write_fasta(out_records, out_fa)
    print(f"\nWrote {len(out_records)} sequences to {out_fa}")

    # --- update SAF (whole-transcript features, same as cds_mane.saf) -------
    out_lines = []
    updated = 0
    with open(args.saf) as f:
        out_lines.append(f.readline())   # header row
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in saf_replacements:
                for new_id, new_len in saf_replacements[parts[0]]:
                    # SAF columns: GeneID, Chr, Start, End, Strand
                    out_lines.append(f"{new_id}\t{new_id}\t{parts[2]}\t{new_len}\t{parts[4]}\n")
                    updated += 1
            else:
                out_lines.append(line)
    if updated != sum(len(v) for v in saf_replacements.values()):
        raise SystemExit("SAF update incomplete -- HLA transcript IDs not found in SAF")
    with open(out_saf, "w") as f:
        f.writelines(out_lines)
    print(f"Wrote updated SAF to {out_saf}  ({updated} HLA row(s))")


if __name__ == "__main__":
    main()
