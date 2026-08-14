#!/usr/bin/env python3
"""
filter_vcfs_to_bed.py

For every VCF found across all tools, platforms, and barcodes, write a
filtered copy retaining only records that overlap the 9-gene BED regions.
Filtered VCFs are written alongside the originals as *_filtered.vcf.
"""

import os
import gzip
import re


# ── Paths ─────────────────────────────────────────────────────────────────────
BASE     = "/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers"
BED_PATH = os.path.join(BASE, "3539131_Covered_DNA_Screen_9genes.bed")

TOOLS     = ["sniffles", "cutesv", "svim", "debreak"]
PLATFORMS = ["ont", "pacbio"]

BARCODE_ORDER = [
    "bc1001", "bc1002", "bc1003", "bc1004",
    "bc1005", "bc1006", "bc1007", "bc1008",
    "bc1097", "bc1098", "bc1099", "bc1100",
    "bc1101", "bc1102", "bc1103", "bc1104",
]


# ── Helpers ───────────────────────────────────────────────────────────────────
def normalise_chrom(c):
    return str(c).strip().lstrip("chr") if c else c


def load_bed(bed_path):
    regions = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            regions.append((normalise_chrom(parts[0]), int(parts[1]), int(parts[2])))
    return regions


def sv_overlaps_bed(sv_chrom, sv_start, sv_end, bed_regions):
    nc = normalise_chrom(sv_chrom)
    for (b_chrom, b_start, b_end) in bed_regions:
        if nc == b_chrom and sv_start <= b_end and sv_end >= b_start:
            return True
    return False


def get_info_tag(info_str, tag):
    m = re.search(rf'(?:^|;){re.escape(tag)}=([^;]+)', info_str)
    return m.group(1) if m else None


def open_vcf_read(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def find_vcf(tool, platform, barcode):
    if tool == "sniffles":
        candidates = [
            os.path.join(BASE, tool, f"{tool}_{platform}", barcode, f"sniffles2_{barcode}.vcf"),
            os.path.join(BASE, tool, f"{tool}_{platform}", barcode, f"sniffles2_{barcode}.vcf.gz"),
        ]
    elif tool == "cutesv":
        candidates = [
            os.path.join(BASE, tool, f"{tool}_{platform}", barcode, f"cutesv_{barcode}.vcf"),
            os.path.join(BASE, tool, f"{tool}_{platform}", barcode, f"cutesv_{barcode}.vcf.gz"),
        ]
    elif tool == "svim":
        candidates = [
            os.path.join(BASE, tool, f"{tool}_{platform}", barcode, "variants.vcf"),
            os.path.join(BASE, tool, f"{tool}_{platform}", barcode, "variants.vcf.gz"),
        ]
    elif tool == "debreak":
        result_dir = os.path.join(BASE, tool, f"{tool}_{platform}", "debreak_results", barcode)
        candidates = []
        if os.path.isdir(result_dir):
            for fname in os.listdir(result_dir):
                if fname.endswith((".vcf", ".vcf.gz")):
                    candidates.append(os.path.join(result_dir, fname))
    else:
        candidates = []

    tool_dir = os.path.join(BASE, tool, f"{tool}_{platform}")
    if os.path.isdir(tool_dir):
        for root, dirs, files in os.walk(tool_dir):
            for fname in files:
                if fname.endswith((".vcf", ".vcf.gz")) and barcode in root:
                    candidates.append(os.path.join(root, fname))

    for path in candidates:
        if os.path.exists(path):
            return path
    return None


def output_path(vcf_path):
    """Derive filtered output path — strips .gz if present, inserts _filtered."""
    base = vcf_path[:-3] if vcf_path.endswith(".gz") else vcf_path
    stem, ext = os.path.splitext(base)
    return f"{stem}_filtered{ext}"


# ── Load BED ──────────────────────────────────────────────────────────────────
bed_regions = load_bed(BED_PATH)
print(f"Loaded {len(bed_regions)} BED regions from {BED_PATH}\n")


# ── Filter loop ───────────────────────────────────────────────────────────────
summary = []

for barcode in BARCODE_ORDER:
    for tool in TOOLS:
        for platform in PLATFORMS:
            vcf_path = find_vcf(tool, platform, barcode)

            if vcf_path is None:
                print(f"  SKIP (no VCF) — {barcode} | {tool} | {platform}")
                continue

            out_path    = output_path(vcf_path)
            total_in    = 0
            total_kept  = 0
            header_lines = []
            data_lines   = []

            try:
                with open_vcf_read(vcf_path) as fh:
                    for line in fh:
                        if line.startswith("#"):
                            header_lines.append(line)
                            continue

                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 8:
                            continue

                        total_in += 1
                        chrom = parts[0]
                        pos   = int(parts[1])
                        info  = parts[7]

                        end_tag = get_info_tag(info, "END")
                        try:
                            sv_end = int(end_tag) if end_tag else pos
                        except ValueError:
                            sv_end = pos

                        if sv_overlaps_bed(chrom, pos, sv_end, bed_regions):
                            total_kept += 1
                            data_lines.append(line)

            except Exception as e:
                print(f"  ERROR reading {vcf_path}: {e}")
                continue

            # Write filtered VCF (plain text regardless of input compression)
            try:
                with open(out_path, "w") as out_fh:
                    out_fh.writelines(header_lines)
                    out_fh.writelines(data_lines)
            except Exception as e:
                print(f"  ERROR writing {out_path}: {e}")
                continue

            print(f"  {barcode} | {tool:8s} | {platform.upper():6s} → "
                  f"{total_kept:4d} / {total_in:6d} records kept  →  {out_path}")

            summary.append({
                "barcode":    barcode,
                "tool":       tool,
                "platform":   platform,
                "total_in":   total_in,
                "kept":       total_kept,
                "out_path":   out_path,
            })

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"Done. {len(summary)} VCFs filtered.")
total_written = sum(r["kept"] for r in summary)
print(f"Total records retained across all filtered VCFs: {total_written}")