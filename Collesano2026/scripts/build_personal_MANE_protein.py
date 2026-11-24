#!/usr/bin/env python3
"""
Build a personalized MANE protein FASTA to match a *_personal_MANE.fa file.

Only HLA-A/B/C differ between a sample's personal_MANE.fa and the reference
MANE mRNA FASTA (build_personal_MANE.py replaces those with sample allele
CDS sequences; since 2026-10-05 build_personal_MANE_utr.py keeps the
MANE UTRs and tags the CDS span as " cds:<start>-<end>" in the header).
Every other transcript is untouched, so
its protein is copied straight from the reference MANE protein FASTA rather
than re-translated (translating full-length UTR-CDS-UTR mRNAs would need
CDS-boundary info that cds_mane.saf does not carry -- it spans whole
transcripts, not just the CDS).

Usage:
    python3 build_personal_MANE_protein.py \
        --personal_mrna HCC1937_personal_MANE.fa \
        --ref_protein   MANE.GRCh38.v1.5.ensembl_protein.faa \
        --out           HCC1937_personal_MANE_protein.faa
"""

import argparse
import re

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}


def parse_fasta(path):
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


def translate_cds(seq):
    """Translate an in-frame CDS (ATG...stop). Returns protein without the stop codon."""
    seq = seq.upper()
    if len(seq) % 3 != 0:
        raise ValueError(f"CDS length {len(seq)} is not a multiple of 3")
    aa = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        residue = CODON_TABLE.get(codon, 'X')
        if residue == '*':
            break
        aa.append(residue)
    return "".join(aa)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--personal_mrna", required=True)
    ap.add_argument("--ref_protein", default="MANE.GRCh38.v1.5.ensembl_protein.faa")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # --- reference protein lookup, keyed by transcript ENST id -------------
    ref_by_enst = {}
    for header, seq in parse_fasta(args.ref_protein):
        m = re.search(r'transcript:(\S+)', header)
        if m:
            ref_by_enst[m.group(1)] = (header, seq)

    out_records = []
    n_copied = 0
    n_translated = 0

    KNOWN_HLA_BASE_IDS = {
        "ENST00000376809.10", "ENST00000412585.7", "ENST00000376228.10",
    }

    for header, seq in parse_fasta(args.personal_mrna):
        transcript_id = header.split()[0]
        # Decide per record: translate this sample's own personalised CDS, or copy the
        # reference protein. Match on the known base HLA transcript ID, NOT on an "allele:"
        # header tag: that tag marks only heterozygous (allele-suffixed) entries, while
        # homozygous genes keep the original ID and header even though the sequence is the
        # sample's own allele. Matching on the tag would copy the generic reference protein
        # for every homozygous personalised gene.
        base_transcript_id = transcript_id if transcript_id in KNOWN_HLA_BASE_IDS \
            else transcript_id.rsplit("_", 1)[0]
        is_allele = base_transcript_id in KNOWN_HLA_BASE_IDS

        if is_allele:
            # sample-specific HLA CDS -- translate directly. HLA records carry the MANE UTRs
            # and a " cds:<start>-<end>" header tag (build_personal_MANE_utr.py); slice the
            # CDS out first.
            # Records without the tag are the old CDS-only format.
            m_cds = re.search(r'\bcds:(\d+)-(\d+)', header)
            if m_cds:
                seq = seq[int(m_cds.group(1)) - 1:int(m_cds.group(2))]
            protein = translate_cds(seq)
            is_suffixed = transcript_id != base_transcript_id   # heterozygous (allele-split) vs homozygous (same ID)

            if is_suffixed and base_transcript_id in ref_by_enst:
                # heterozygous: construct the allele-suffixed protein ID, same as before
                ref_header = ref_by_enst[base_transcript_id][0]
                ref_ensp = ref_header.split()[0]
                allele_suffix = transcript_id.rsplit("_", 1)[1]
                ensp_id = f"{ref_ensp}_{allele_suffix}"
                rest = ref_header[len(ref_ensp):]
                rest = re.sub(r'transcript:\S+', f'transcript:{transcript_id}', rest)
            elif (not is_suffixed) and transcript_id in ref_by_enst:
                # homozygous: SAME (unsuffixed) protein ID as the reference/copy path
                # would have produced -- downstream files already reference this ID.
                ref_header = ref_by_enst[transcript_id][0]
                ref_ensp = ref_header.split()[0]
                ensp_id = ref_ensp
                rest = ref_header[len(ref_ensp):]
            else:
                ensp_id = transcript_id
                rest = header[len(transcript_id):]
            new_header = f"{ensp_id}{rest}"
            out_records.append((new_header, protein))
            n_translated += 1
            print(f"  translated {transcript_id}: {len(seq)} nt -> {len(protein)} aa"
                  + ("" if is_suffixed else "  [homozygous, ID unchanged]"))
        else:
            ref = ref_by_enst.get(transcript_id)
            if ref is None:
                print(f"  WARNING: no reference protein found for {transcript_id} (non-coding / no ENSP?) -- skipping")
                continue
            out_records.append(ref)
            n_copied += 1

    write_fasta(out_records, args.out)
    print(f"\nWrote {len(out_records)} protein sequences to {args.out}")
    print(f"  {n_copied} copied from reference, {n_translated} translated from sample alleles")


if __name__ == "__main__":
    main()
