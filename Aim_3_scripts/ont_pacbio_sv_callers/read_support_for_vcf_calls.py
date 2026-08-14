#!/usr/bin/env python3

import os
import gzip
import re
import pandas as pd
from statistics import mean, median


# ── Paths ─────────────────────────────────────────────────────────────────────
BASE                  = "/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers"
OUTPUT_TSV            = os.path.join(BASE, "sv_read_support_distribution.tsv")
REAL_BREAKPOINTS_TSV  = os.path.join(BASE, "real_breakpoints.tsv")
BED_PATH              = os.path.join(BASE, "3539131_Covered_DNA_Screen_9genes.bed")

TOOLS     = ["sniffles", "cutesv", "svim", "debreak"]
PLATFORMS = ["ont", "pacbio"]

TOOL_LABELS = {
    "sniffles": "Sniffles (default run)",
    "cutesv":   "CuteSV (run with suggested params)",
    "svim":     "SVIM (default run)",
    "debreak":  "DeBreak (default run)",
}

BARCODE_ORDER = [
    "bc1001", "bc1002", "bc1003", "bc1004",
    "bc1005", "bc1006", "bc1007", "bc1008",
    "bc1097", "bc1098", "bc1099", "bc1100",
    "bc1101", "bc1102", "bc1103", "bc1104",
]

APPROX_BARCODES = {
    "bc1006",   # PCSK9 gene window (GRCh38: chr1:55,039,476-55,064,853)
    "bc1103",   # MLH1  gene window (GRCh38: chr3:36,993,226-37,050,896)
}
APPROX_GENE_LABELS = {
    "bc1006": "PCSK9",
    "bc1103": "MLH1",
}

TYPE_MAP = {"Deletion": "DEL", "Duplication": "DUP"}


# ── Helpers ───────────────────────────────────────────────────────────────────
def normalise_chrom(c):
    return str(c).strip().lstrip("chr") if c else c


def parse_coord(coord_str):
    s = str(coord_str).strip()
    if s in ("NA", "", "nan") or pd.isna(coord_str):
        return None, None
    try:
        s = s.replace(",", "")
        chrom, pos_str = s.split(":", 1)
        return normalise_chrom(chrom.strip()), int(pos_str.strip())
    except (ValueError, AttributeError):
        return None, None


def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def get_info_tag(info_str, tag):
    m = re.search(rf'(?:^|;){re.escape(tag)}=([^;]+)', info_str)
    return m.group(1) if m else None


def get_format_val(fmt_str, samp_str, tag):
    if not fmt_str or not samp_str:
        return None
    d = dict(zip(fmt_str.split(":"), samp_str.split(":")))
    return d.get(tag)


def parse_read_support(info, fmt, samp, tool):
    try:
        if tool == "sniffles":
            val = get_info_tag(info, "SUPPORT") or get_format_val(fmt, samp, "DV")
        elif tool == "cutesv":
            val = get_info_tag(info, "RE") or get_format_val(fmt, samp, "DV")
        elif tool == "svim":
            val = get_info_tag(info, "SUPPORT")
            if val is None:
                ad = get_format_val(fmt, samp, "AD")
                if ad and "," in ad:
                    val = ad.split(",")[1]
        elif tool == "debreak":
            val = get_info_tag(info, "SUPPREAD")
        else:
            val = None
        return int(val) if val is not None else None
    except (ValueError, TypeError):
        return None


def sv_overlaps_window(sv_chrom, sv_start, sv_end, win_chrom, win_start, win_end):
    return (normalise_chrom(sv_chrom) == normalise_chrom(win_chrom)
            and sv_start <= win_end
            and sv_end   >= win_start)


# ── BED loading and filter ────────────────────────────────────────────────────
def load_bed(bed_path):
    regions = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
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


# ── Load real breakpoints ─────────────────────────────────────────────────────
real_bp           = {}
barcode_to_sample = {}

if os.path.exists(REAL_BREAKPOINTS_TSV):
    bp_df = pd.read_csv(REAL_BREAKPOINTS_TSV, sep="\t", dtype=str)
    bp_df.columns = bp_df.columns.str.strip()
    for _, bprow in bp_df.iterrows():
        raw_bc    = str(bprow.get("Barcode", "")).strip()
        bc_num    = raw_bc.lstrip("bc").zfill(4)
        barcode   = f"bc{bc_num}"
        sample_id = str(bprow.get("Sample ID", "NA")).strip()
        cnv_type  = str(bprow.get("CNV event based on manual IGV inspection",
                                   bprow.get("CNV event based on LRS", ""))).strip()
        ts_chrom, ts_pos = parse_coord(bprow.get("True CNV Start", "NA"))
        te_chrom, te_pos = parse_coord(bprow.get("True CNV End",   "NA"))
        real_bp[barcode] = {
            "sample_id":        sample_id,
            "cnv_type":         cnv_type,
            "expected_svtype":  TYPE_MAP.get(cnv_type),
            "true_start_chrom": ts_chrom,
            "true_start_pos":   ts_pos,
            "true_end_chrom":   te_chrom,
            "true_end_pos":     te_pos,
        }
        barcode_to_sample[barcode] = sample_id
    print(f"Loaded real breakpoints for {len(real_bp)} barcodes.")
else:
    print(f"WARNING: real_breakpoints.tsv not found at {REAL_BREAKPOINTS_TSV}")

bed_regions = load_bed(BED_PATH)
print(f"Loaded {len(bed_regions)} BED regions from {BED_PATH}")


# ── VCF path finder ───────────────────────────────────────────────────────────
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


# ── CNV concordance assessment ────────────────────────────────────────────────
def assess_cnv_concordance(svs, barcode, bp_info):
    cnv_type        = bp_info["cnv_type"]
    expected_svtype = bp_info["expected_svtype"]
    true_chrom      = bp_info["true_start_chrom"]
    true_start      = bp_info["true_start_pos"]
    true_end        = bp_info["true_end_pos"]
    is_approx       = barcode in APPROX_BARCODES
    note_suffix     = (f" (approx gene-window: {APPROX_GENE_LABELS[barcode]})"
                       if is_approx else "")

    no_cnv_row = {
        "cnv_type":               cnv_type,
        "true_start":             "NA",
        "true_end":               "NA",
        "best_match_chrom":       "NA",
        "best_match_start":       "NA",
        "best_match_end":         "NA",
        "best_match_svtype":      "NA",
        "best_match_read_support":"NA",
        "start_offset":           "NA",
        "end_offset":             "NA",
        "concordance_note":       "No true CNV / NA barcode",
    }

    if expected_svtype is None or true_start is None:
        return no_cnv_row

    if is_approx:
        candidates = [sv for sv in svs
                      if sv["svtype"] == expected_svtype
                      and sv_overlaps_window(sv["chrom"], sv["start"], sv["end"],
                                             true_chrom, true_start, true_end)]
        if not candidates:
            candidates = [sv for sv in svs
                          if sv["svtype"] == expected_svtype
                          and normalise_chrom(sv["chrom"]) == normalise_chrom(true_chrom)]
    else:
        candidates = [sv for sv in svs
                      if sv["svtype"] == expected_svtype
                      and normalise_chrom(sv["chrom"]) == normalise_chrom(true_chrom)]

    if not candidates:
        return {
            "cnv_type":               cnv_type,
            "true_start":             true_start,
            "true_end":               true_end,
            "best_match_chrom":       true_chrom,
            "best_match_start":       "NA",
            "best_match_end":         "NA",
            "best_match_svtype":      "NA",
            "best_match_read_support":"NA",
            "start_offset":           "NA",
            "end_offset":             "NA",
            "concordance_note":       (f"No {expected_svtype} found on "
                                       f"chr{true_chrom} in VCF{note_suffix}"),
        }

    if is_approx:
        best = max(candidates, key=lambda sv:
                   min(sv["end"], true_end) - max(sv["start"], true_start))
    else:
        best = min(candidates, key=lambda sv:
                   abs(sv["start"] - true_start) + abs(sv["end"] - true_end))

    return {
        "cnv_type":               cnv_type,
        "true_start":             true_start,
        "true_end":               true_end,
        "best_match_chrom":       best["chrom"],
        "best_match_start":       best["start"],
        "best_match_end":         best["end"],
        "best_match_svtype":      best["svtype"],
        "best_match_read_support":best.get("support", "NA"),
        "start_offset":           (best["start"] - true_start
                                   if isinstance(best["start"], int) else "NA"),
        "end_offset":             (best["end"] - true_end
                                   if isinstance(best["end"], int) else "NA"),
        "concordance_note":       f"OK{note_suffix}",
    }


# ── Main loop ─────────────────────────────────────────────────────────────────
NA_VALS = {
    "Total SV calls":            "N/A",
    "Min read support":          "N/A",
    "Max read support":          "N/A",
    "Mean read support":         "N/A",
    "Median read support":       "N/A",
    "Min SV size (bp)":          "N/A",
    "Max SV size (bp)":          "N/A",
    "Median SV size (bp)":       "N/A",
    "CNV type (truth)":          "N/A",
    "True CNV start":            "N/A",
    "True CNV end":              "N/A",
    "Best match chrom":          "N/A",
    "Best match start":          "N/A",
    "Best match end":            "N/A",
    "Best match SVTYPE":         "N/A",
    "Best match read support":   "N/A",
    "Start offset (bp)":         "N/A",
    "End offset (bp)":           "N/A",
    "Concordance note":          "VCF missing",
    "Max read support == overlapping call?": "N/A",
}

rows = []

for barcode in BARCODE_ORDER:
    bp_info   = real_bp.get(barcode)
    sample_id = bp_info["sample_id"] if bp_info else "NA"
    has_real_bp = (
        bp_info is not None
        and bp_info["true_start_pos"] is not None
        and bp_info["true_end_pos"] is not None
    )

    for tool in TOOLS:
        tool_label = TOOL_LABELS[tool]
        for platform in PLATFORMS:
            plat_label = platform.upper()
            vcf_path   = find_vcf(tool, platform, barcode)

            row = {
                "Barcode":   barcode,
                "Sample ID": sample_id,
                "Tool":      tool_label,
                "Platform":  plat_label,
            }

            if vcf_path is None:
                row.update(NA_VALS)
                rows.append(row)
                continue

            # ── Parse VCF, keeping only BED-overlapping SVs ───────────────────
            support_vals = []
            size_vals    = []
            svs          = []

            try:
                with open_vcf(vcf_path) as fh:
                    for line in fh:
                        if line.startswith("#"):
                            continue
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 8:
                            continue

                        chrom = parts[0]
                        pos   = int(parts[1])
                        info  = parts[7]
                        fmt   = parts[8] if len(parts) > 8 else ""
                        samp  = parts[9] if len(parts) > 9 else ""

                        end_tag = get_info_tag(info, "END")
                        try:
                            sv_end = int(end_tag) if end_tag else pos
                        except ValueError:
                            sv_end = pos

                        if not sv_overlaps_bed(chrom, pos, sv_end, bed_regions):
                            continue

                        svtype = (get_info_tag(info, "SVTYPE") or ".").upper()
                        rs     = parse_read_support(info, fmt, samp, tool)
                        sv_size = sv_end - pos  # bp size of the SV call

                        if rs is not None:
                            support_vals.append(rs)

                        if sv_size > 0:
                            size_vals.append(sv_size)

                        svs.append({
                            "svtype":  svtype,
                            "chrom":   normalise_chrom(chrom),
                            "start":   pos,
                            "end":     sv_end,
                            "support": rs,
                            "size":    sv_size,
                        })

            except Exception as e:
                print(f"  ERROR reading {vcf_path}: {e}")
                row.update({k: f"ERROR: {e}" if k != "Concordance note" else "ERROR"
                             for k in NA_VALS})
                rows.append(row)
                continue

            # ── Read support summary ──────────────────────────────────────────
            total_calls = len(svs)
            if support_vals:
                max_support = max(support_vals)
                row["Total SV calls"]      = total_calls
                row["Min read support"]    = min(support_vals)
                row["Max read support"]    = max_support
                row["Mean read support"]   = round(mean(support_vals), 2)
                row["Median read support"] = median(support_vals)
            else:
                max_support = None
                row["Total SV calls"]      = total_calls
                row["Min read support"]    = "N/A"
                row["Max read support"]    = "N/A"
                row["Mean read support"]   = "N/A"
                row["Median read support"] = "N/A"

            # ── SV size summary ───────────────────────────────────────────────
            if size_vals:
                row["Min SV size (bp)"]    = min(size_vals)
                row["Max SV size (bp)"]    = max(size_vals)
                row["Median SV size (bp)"] = median(size_vals)
            else:
                row["Min SV size (bp)"]    = "N/A"
                row["Max SV size (bp)"]    = "N/A"
                row["Median SV size (bp)"] = "N/A"

            # ── CNV concordance ───────────────────────────────────────────────
            if not has_real_bp:
                conc = {k: "No real BP defined" for k in [
                    "cnv_type", "true_start", "true_end",
                    "best_match_chrom", "best_match_start", "best_match_end",
                    "best_match_svtype", "best_match_read_support",
                    "start_offset", "end_offset", "concordance_note",
                ]}
                max_eq_overlap = "N/A"
            else:
                conc = assess_cnv_concordance(svs, barcode, bp_info)

                bm_rs = conc["best_match_read_support"]
                if (max_support is not None
                        and isinstance(bm_rs, int)
                        and conc["concordance_note"].startswith("OK")):
                    max_eq_overlap = "YES" if max_support == bm_rs else "NO"
                elif conc["concordance_note"].startswith("No "):
                    max_eq_overlap = "NO"
                else:
                    max_eq_overlap = "N/A"

            row["CNV type (truth)"]          = conc["cnv_type"]
            row["True CNV start"]            = conc["true_start"]
            row["True CNV end"]              = conc["true_end"]
            row["Best match chrom"]          = conc["best_match_chrom"]
            row["Best match start"]          = conc["best_match_start"]
            row["Best match end"]            = conc["best_match_end"]
            row["Best match SVTYPE"]         = conc["best_match_svtype"]
            row["Best match read support"]   = conc["best_match_read_support"]
            row["Start offset (bp)"]         = conc["start_offset"]
            row["End offset (bp)"]           = conc["end_offset"]
            row["Concordance note"]          = conc["concordance_note"]
            row["Max read support == overlapping call?"] = max_eq_overlap

            rows.append(row)
            print(f"  {barcode} | {tool_label} | {plat_label} → "
                  f"{total_calls} BED-overlapping calls | {conc['concordance_note']} | "
                  f"max==best: {max_eq_overlap}")


# ── Write output ──────────────────────────────────────────────────────────────
COLUMNS = [
    "Barcode", "Sample ID", "Tool", "Platform",
    "Total SV calls",
    "Min read support", "Max read support",
    "Mean read support", "Median read support",
    "Min SV size (bp)", "Max SV size (bp)", "Median SV size (bp)",
    "CNV type (truth)",
    "True CNV start", "True CNV end",
    "Best match chrom", "Best match start", "Best match end",
    "Best match SVTYPE", "Best match read support",
    "Start offset (bp)", "End offset (bp)",
    "Concordance note",
    "Max read support == overlapping call?",
]

out_df = pd.DataFrame(rows, columns=COLUMNS)
out_df.to_csv(OUTPUT_TSV, sep="\t", index=False)
print(f"\nWritten : {OUTPUT_TSV}")
print(f"Shape   : {out_df.shape[0]} rows x {out_df.shape[1]} columns")