#!/usr/bin/env python3

import os
import argparse
import pandas as pd

# ── Mode argument ─────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Compare DeBreak CNV calls to truth set.")
parser.add_argument("--mode", choices=["ont", "ont_round2", "pacbio"], required=True,
                    help="Sequencing mode: 'ont' or 'pacbio'")
args = parser.parse_args()

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR = "/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/debreak"
TSV_FILE = "/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/real_breakpoints.tsv"

if args.mode == "ont":
    VCF_BASE   = os.path.join(BASE_DIR, "debreak_ont")
    OUTPUT_TSV = os.path.join(VCF_BASE, "debreak_cnv_comparison.tsv")
elif args.mode == "ont_round2":
    VCF_BASE   = os.path.join(BASE_DIR, "debreak_ont_round2")
    OUTPUT_TSV = os.path.join(VCF_BASE, "debreak_cnv_comparison_round2.tsv")
elif args.mode == "pacbio":
    # NOTE: assumed to mirror the ONT naming (debreak_pacbio); confirm this
    # matches your actual PacBio DeBreak output folder if it differs.
    VCF_BASE   = os.path.join(BASE_DIR, "debreak_pacbio")
    OUTPUT_TSV = os.path.join(VCF_BASE, "debreak_cnv_comparison.tsv")

print(f"Mode     : {args.mode}")
print(f"VCF base : {VCF_BASE}")
print(f"Output   : {OUTPUT_TSV}")

# ── Barcodes with approximate/gene-window coordinates ────────────────────────
APPROX_BARCODES = {
    "bc1006",   # PCSK9 gene window (GRCh38: chr1:55,035,547-55,065,852)
    "bc1103",   # MLH1  gene window (GRCh38: chr3:36,989,517-37,051,847)
    "bc1116",   # BRCA2 gene window (GRCh38: chr13:32,311,507-32,401,269)
    "bc1117",   # BRCA1 gene window (GRCh38: chr17:43,043,294-43,129,365)
    "bc1118",   # PCSK9 gene window (GRCh38: chr1:55,035,547-55,065,852)
    "bc1119",   # MSH2  gene window (GRCh38: chr2:47,399,155-47,484,176)
    "bc1125",   # BRCA2 gene window (GRCh38: chr13:32,311,507-32,401,269)
    "bc1126",   # BRCA2 gene window (GRCh38: chr13:32,311,507-32,401,269)
}

APPROX_GENE_LABELS = {
    "bc1006": "PCSK9",
    "bc1103": "MLH1",
    "bc1116": "BRCA2",
    "bc1117": "BRCA1",
    "bc1118": "PCSK9",
    "bc1119": "MSH2",
    "bc1125": "BRCA2",
    "bc1126": "BRCA2",
}

# ── Load truth TSV ────────────────────────────────────────────────────────────
truth = pd.read_csv(TSV_FILE, sep="\t")
truth.columns = [c.strip() for c in truth.columns]

truth = truth.rename(columns={
    "Barcode": "barcode",
    "Sample ID": "sample_id",
    "CNV event based on LRS": "cnv_type",
    "True CNV Start": "true_start_raw",
    "True CNV End": "true_end_raw"
})


def normalise_chrom(chrom):
    """Strip 'chr' prefix for consistent comparison."""
    return str(chrom).strip().lstrip("chr") if chrom else chrom


def parse_coord(coord_str):
    """Parse 'chr:pos' string -> (chrom, int pos). Returns (None, None) for NA."""
    if pd.isna(coord_str) or str(coord_str).strip().upper() == "NA":
        return None, None
    coord_str = str(coord_str).replace(",", "").strip()
    chrom, pos = coord_str.split(":")
    return normalise_chrom(chrom), int(pos.strip())


def parse_vcf(vcf_path):
    """Parse a DeBreak VCF and return list of dicts with type, chrom, start, end.

    DeBreak only emits an explicit END= tag for symbolic ALT records (e.g.
    <DUP>). DEL/INS records with sequence-based REF/ALT instead carry
    SVLEN=, so END must be derived as POS + abs(SVLEN) when END is absent.
    """
    svs = []
    if not os.path.exists(vcf_path):
        print(f"  WARNING: VCF not found: {vcf_path}")
        return svs
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom = normalise_chrom(parts[0])
            pos   = int(parts[1])
            info  = parts[7]

            svtype, end, svlen = None, None, None
            for field in info.split(";"):
                if field.startswith("SVTYPE="):
                    svtype = field.split("=")[1].upper()
                elif field.startswith("END="):
                    try:
                        end = int(field.split("=")[1])
                    except ValueError:
                        pass
                elif field.startswith("SVLEN="):
                    try:
                        svlen = int(field.split("=")[1])
                    except ValueError:
                        pass

            if svtype is None:
                continue

            if end is None and svlen is not None:
                end = pos + abs(svlen)

            if end is None:
                continue

            svs.append({"svtype": svtype, "chrom": chrom, "start": pos, "end": end})
    return svs


def overlaps(sv, chrom, window_start, window_end):
    """Return True if the SV overlaps the given genomic window."""
    return sv["chrom"] == chrom and sv["start"] <= window_end and sv["end"] >= window_start


# ── Map CNV type strings to VCF SVTYPE ───────────────────────────────────────
TYPE_MAP = {"Deletion": "DEL", "Duplication": "DUP"}

results = []

for _, row in truth.iterrows():
    bc       = str(row["barcode"]).zfill(4)
    barcode  = f"bc{bc}"
    cnv_type = str(row["cnv_type"]).strip()
    expected_svtype = TYPE_MAP.get(cnv_type)
    is_approx   = barcode in APPROX_BARCODES
    note_suffix = f" (approx gene-window: {APPROX_GENE_LABELS[barcode]})" if is_approx else ""

    true_chrom_s, true_start = parse_coord(row["true_start_raw"])
    true_chrom_e, true_end   = parse_coord(row["true_end_raw"])

    if expected_svtype is None or true_start is None:
        results.append({
            "barcode": barcode,
            "sample_id": row["sample_id"],
            "cnv_type": cnv_type,
            "true_start": row["true_start_raw"],
            "true_end": row["true_end_raw"],
            "best_match_chrom": "NA",
            "best_match_start": "NA",
            "best_match_end": "NA",
            "start_offset": "NA",
            "end_offset": "NA",
            "note": "No true CNV / NA barcode"
        })
        continue

    # DeBreak output: <VCF_BASE>/<barcode>/debreak.vcf (fixed filename)
    vcf_path = os.path.join(VCF_BASE, barcode, "debreak.vcf")
    svs      = parse_vcf(vcf_path)

    if is_approx:
        candidates = [
            sv for sv in svs
            if sv["svtype"] == expected_svtype
            and overlaps(sv, true_chrom_s, true_start, true_end)
        ]
        if not candidates:
            candidates = [
                sv for sv in svs
                if sv["svtype"] == expected_svtype and sv["chrom"] == true_chrom_s
            ]
    else:
        candidates = [
            sv for sv in svs
            if sv["svtype"] == expected_svtype and sv["chrom"] == true_chrom_s
        ]

    if not candidates:
        results.append({
            "barcode": barcode,
            "sample_id": row["sample_id"],
            "cnv_type": cnv_type,
            "true_start": true_start,
            "true_end": true_end,
            "best_match_chrom": true_chrom_s,
            "best_match_start": "NA",
            "best_match_end": "NA",
            "start_offset": "NA",
            "end_offset": "NA",
            "note": f"No {expected_svtype} found on chr{true_chrom_s} in VCF{note_suffix}"
        })
        continue

    if is_approx:
        best = max(candidates, key=lambda sv:
                   min(sv["end"], true_end) - max(sv["start"], true_start))
    else:
        best = min(candidates, key=lambda sv:
                   abs(sv["start"] - true_start) + abs(sv["end"] - true_end))

    results.append({
        "barcode": barcode,
        "sample_id": row["sample_id"],
        "cnv_type": cnv_type,
        "true_start": true_start,
        "true_end": true_end,
        "best_match_chrom": best["chrom"],
        "best_match_start": best["start"],
        "best_match_end": best["end"],
        "start_offset": best["start"] - true_start,
        "end_offset": best["end"]   - true_end,
        "note": f"OK{note_suffix}"
    })

# ── Write output ──────────────────────────────────────────────────────────────
out_df = pd.DataFrame(results)
out_df.to_csv(OUTPUT_TSV, sep="\t", index=False)
print(f"\nWritten: {OUTPUT_TSV}")
print(out_df.to_string(index=False))
