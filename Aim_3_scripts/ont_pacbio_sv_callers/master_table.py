#!/usr/bin/env python3

import os
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = "/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers"
OUTPUT_TSV = os.path.join(BASE, "master_cnv_comparison.tsv")

SOFTCLIP_RESULTS_DIR = (
    "/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/"
    "ONT_PacBio_CNV_calling/automated_precise_breakpoint_finding/"
    "batch_scripts/results/sv_analysis_py_all_bc"
)

REAL_BREAKPOINTS_TSV = os.path.join(BASE, "real_breakpoints.tsv")

TOOLS = ["sniffles", "cutesv", "svim", "debreak"]
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

# ── aCGH reference coordinates (GRCh38, chrom without 'chr') ─────────────────
ACGH_COORDS = {
    "DNS000601": ("13", 32325998, 32326724),
    "DNS002261": ("3",  37046499, 37050922),
    "DNS004024": ("17", 43111887, 43117024),
    "DNS004144": ("2",  20900554, 21411082),
    "DNS004937": ("17", 43045165, 43045883),
    "DNS006443": ("16", 23606156, 23614296),
    "DNS006683": ("1",  55037565, 55040505),
    "DNS006916": ("17", 43078723, 43091151),
    "DNS008387": ("19", 11096325, 11166174),
    "DNS010349": ("2",  20997587, 21027015),
    # Inconclusive aCGH — coordinates exist but result is flagged as inconclusive
    "DNS004681": ("1",  55040596, 55060215),
    "DNS007821": ("3",  36993256, 37012096),
}

# Samples whose concordance column should always read "Inconclusive aCGH data"
INCONCLUSIVE_ACGH = {"DNS004681", "DNS007821"}

# For these samples, recompute offsets against the aCGH coords rather than the
# gene-window true_start/true_end stored in the comparison TSV.
OFFSET_OVERRIDE = {
    "DNS004681": ("1", 55040596, 55060215),
    "DNS007821": ("3", 36993256, 37012096),
}

# ── Load real breakpoints TSV ─────────────────────────────────────────────────
# Columns: Barcode | Sample ID | CNV event based on LRS | True CNV Start | True CNV End
# True CNV Start/End format: "chrom:pos,pos,pos" (commas as thousand separators)
def parse_bp_coord(coord_str):
    """Parse 'chrom:pos' where pos may contain comma thousand-separators.
    Returns (chrom_str, int_pos) or (None, None) on failure / NA.
    """
    s = str(coord_str).strip()
    if s in ("NA", "", "nan"):
        return None, None
    try:
        chrom, pos_str = s.split(":", 1)
        pos = int(pos_str.replace(",", ""))
        return chrom.strip(), pos
    except (ValueError, AttributeError):
        return None, None

real_bp = {}  # barcode_num (str) -> {"cnv_event": ..., "true_start": (chrom, pos), "true_end": (chrom, pos)}
if os.path.exists(REAL_BREAKPOINTS_TSV):
    bp_df = pd.read_csv(REAL_BREAKPOINTS_TSV, sep="\t", dtype=str)
    bp_df.columns = bp_df.columns.str.strip()
    for _, bprow in bp_df.iterrows():
        bc_num = str(bprow.get("Barcode", "")).strip()
        cnv_event = str(bprow.get("CNV event based on LRS", "")).strip()
        ts_chrom, ts_pos = parse_bp_coord(bprow.get("True CNV Start", "NA"))
        te_chrom, te_pos = parse_bp_coord(bprow.get("True CNV End", "NA"))
        real_bp[bc_num] = {
            "cnv_event": cnv_event,
            "true_start_chrom": ts_chrom,
            "true_start_pos":   ts_pos,
            "true_end_chrom":   te_chrom,
            "true_end_pos":     te_pos,
        }
    print(f"Loaded real breakpoints for {len(real_bp)} barcodes.")
else:
    print(f"WARNING: real_breakpoints.tsv not found at {REAL_BREAKPOINTS_TSV}")

# ── Load soft-clip results (batch_mh_results.csv) ─────────────────────────────
# Folder pattern: {barcode_num}_{PLATFORM_CAPITALISED}  e.g. 1001_ONT, 1001_PacBio
# We index by (barcode_num, platform_lower) -> dict of first row values
softclip_data = {}  # (bc_num, platform_lower) -> {"left_clip_pos": val, "right_clip_pos": val}

PLATFORM_FOLDER_LABELS = {
    "ont":    ["ONT"],
    "pacbio": ["PacBio"],
}

for bc_full in BARCODE_ORDER:
    bc_num = bc_full.replace("bc", "")  # "bc1001" -> "1001"
    for plat_lower, folder_suffixes in PLATFORM_FOLDER_LABELS.items():
        for suffix in folder_suffixes:
            folder_name = f"{bc_num}_{suffix}"
            csv_path = os.path.join(SOFTCLIP_RESULTS_DIR, folder_name, "batch_mh_results.csv")
            if os.path.exists(csv_path):
                try:
                    sc_df = pd.read_csv(csv_path, dtype=str)
                    sc_df.columns = sc_df.columns.str.strip()
                    if not sc_df.empty:
                        softclip_data[(bc_num, plat_lower)] = sc_df.iloc[0].to_dict()
                    else:
                        softclip_data[(bc_num, plat_lower)] = {}
                except Exception as e:
                    print(f"  WARNING: Could not read {csv_path}: {e}")
                    softclip_data[(bc_num, plat_lower)] = {}
                break  # found for this platform

# ── Helpers ───────────────────────────────────────────────────────────────────
def is_na_val(val):
    if val is None:
        return True
    if isinstance(val, float) and pd.isna(val):
        return True
    return str(val).strip() in ("NA", "nan", "")

def norm_chrom(chrom):
    return str(chrom).strip().lstrip("chr")

def format_coords(row):
    chrom = row.get("best_match_chrom", "NA")
    start = row.get("best_match_start", "NA")
    end   = row.get("best_match_end",   "NA")
    if is_na_val(start):
        return "did not call"
    try:
        return f"{chrom}:{int(float(start))}-{int(float(end))}"
    except (ValueError, TypeError):
        return "did not call"

def compute_offsets(row, sample_id, note, is_na_barcode):
    """
    Return (start_offset, end_offset) as formatted values.
    For samples in OFFSET_OVERRIDE, recompute against aCGH coords instead of
    the gene-window true_start/true_end from the TSV.
    """
    if is_na_barcode:
        return "N/A", "N/A"

    start = row.get("best_match_start", "NA")
    end   = row.get("best_match_end",   "NA")

    if is_na_val(start):
        placeholder = (
            "Suspected but did not find anything"
            if "approx gene-window" in str(note) else "NO"
        )
        return placeholder, placeholder

    if sample_id in OFFSET_OVERRIDE:
        _, ref_start, ref_end = OFFSET_OVERRIDE[sample_id]
        try:
            s_off = int(float(start)) - ref_start
            e_off = int(float(end))   - ref_end
            return s_off, e_off
        except (ValueError, TypeError):
            return "NO", "NO"

    def _fmt(val):
        if is_na_val(val):
            return "Suspected but did not find anything" if "approx gene-window" in str(note) else "NO"
        try:
            return int(float(val))
        except (ValueError, TypeError):
            return val

    return _fmt(row.get("start_offset")), _fmt(row.get("end_offset"))

def get_concordant(row, sample_id):
    note  = str(row.get("note", ""))
    start = row.get("best_match_start", "NA")

    if "No true CNV / NA barcode" in note:
        return "N/A"

    if sample_id in INCONCLUSIVE_ACGH:
        return "Inconclusive aCGH data"

    if sample_id not in ACGH_COORDS:
        return "No aCGH data"

    acgh_chrom, acgh_start, acgh_end = ACGH_COORDS[sample_id]

    if is_na_val(start):
        if "approx gene-window" in note:
            return "Suspected but did not find anything"
        return "NO"

    try:
        call_chrom = norm_chrom(row.get("best_match_chrom", "NA"))
        call_start = int(float(start))
        call_end   = int(float(row.get("best_match_end", "NA")))
    except (ValueError, TypeError):
        return "NO"

    if call_chrom != acgh_chrom:
        return "NO"

    return "YES" if (call_start <= acgh_end and call_end >= acgh_start) else "NO"

def compute_softclip_offsets(bc_num, plat_lower, is_na_barcode):
    """
    Return (left_clip_offset, right_clip_offset) of soft-clip positions
    against the real breakpoints for this barcode.
    Only computed for true CNVs (not 'No true CNV / NA barcode').
    Returns ("N/A", "N/A") for NA barcodes or missing data.
    """
    if is_na_barcode:
        return "N/A", "N/A"

    bp_info = real_bp.get(bc_num)
    if bp_info is None:
        return "No real BP data", "No real BP data"

    # Skip if no true CNV breakpoints are defined
    if bp_info["true_start_pos"] is None or bp_info["true_end_pos"] is None:
        return "N/A", "N/A"

    sc_row = softclip_data.get((bc_num, plat_lower))
    if sc_row is None:
        return "File missing", "File missing"
    if not sc_row:
        return "Empty file", "Empty file"

    left_raw  = sc_row.get("left_clip_pos",  None)
    right_raw = sc_row.get("right_clip_pos", None)

    def clip_offset(clip_val, ref_pos):
        if clip_val is None or is_na_val(clip_val):
            return "Not found"
        try:
            return int(float(str(clip_val).replace(",", ""))) - ref_pos
        except (ValueError, TypeError):
            return "Parse error"

    left_off  = clip_offset(left_raw,  bp_info["true_start_pos"])
    right_off = clip_offset(right_raw, bp_info["true_end_pos"])
    return left_off, right_off

# ── Load all tool TSVs ────────────────────────────────────────────────────────
data = {}
for tool in TOOLS:
    for platform in PLATFORMS:
        path = os.path.join(BASE, tool, f"{tool}_{platform}", f"{tool}_cnv_comparison.tsv")
        if os.path.exists(path):
            df = pd.read_csv(path, sep="\t", dtype=str)
            df["barcode"] = df["barcode"].astype(str).str.strip()
            data[(tool, platform)] = df.set_index("barcode")
            print(f"  Loaded : {path}")
        else:
            print(f"  MISSING: {path}")
            data[(tool, platform)] = None

# ── Get base info from first available file ───────────────────────────────────
def get_base_info(bc):
    for tool in TOOLS:
        for platform in PLATFORMS:
            df = data.get((tool, platform))
            if df is not None and bc in df.index:
                r = df.loc[bc].to_dict()
                return {
                    "Barcode":    bc,
                    "Sample ID":  r.get("sample_id",  "NA"),
                    "CNV Type":   r.get("cnv_type",   "NA"),
                    "True Start": r.get("true_start",  "NA"),
                    "True End":   r.get("true_end",    "NA"),
                }
    return {"Barcode": bc, "Sample ID": "NA", "CNV Type": "NA",
            "True Start": "NA", "True End": "NA"}

# ── Build master table ────────────────────────────────────────────────────────
rows = []

for bc in BARCODE_ORDER:
    base      = get_base_info(bc)
    row       = dict(base)
    sample_id = base["Sample ID"]
    bc_num    = bc.replace("bc", "")   # numeric portion, e.g. "1001"

    for tool in TOOLS:
        label = TOOL_LABELS[tool]
        for platform in PLATFORMS:
            plat_label = platform.upper()
            df = data.get((tool, platform))

            col_coords   = f"{label} Coordinates on {plat_label}"
            col_start    = f"{label} Offset of Start coord vs True start ({plat_label})"
            col_end      = f"{label} Offset of End coord vs True End ({plat_label})"
            col_concord  = f"{label} on {plat_label} concordant with aCGH?"

            if df is None or bc not in df.index:
                row[col_coords]  = "N/A (file missing)"
                row[col_start]   = "N/A"
                row[col_end]     = "N/A"
                row[col_concord] = "N/A"
                continue

            r       = df.loc[bc].to_dict()
            note    = str(r.get("note", ""))
            is_na_bc = "No true CNV / NA barcode" in note

            s_off, e_off = compute_offsets(r, sample_id, note, is_na_bc)

            row[col_coords]  = "N/A" if is_na_bc else format_coords(r)
            row[col_start]   = s_off
            row[col_end]     = e_off
            row[col_concord] = get_concordant(r, sample_id)

    # ── Soft-clip offset columns (one pair per platform) ──────────────────────
    # Determine is_na_bc from any available tool for this barcode
    # (note field is consistent across tools for the same barcode)
    def _get_note_for_bc(bc_key):
        for t in TOOLS:
            for p in PLATFORMS:
                df2 = data.get((t, p))
                if df2 is not None and bc_key in df2.index:
                    return str(df2.loc[bc_key].to_dict().get("note", ""))
        return ""

    bc_note   = _get_note_for_bc(bc)
    is_na_bc  = "No true CNV / NA barcode" in bc_note

    for platform in PLATFORMS:
        plat_label = platform.upper()
        left_off, right_off = compute_softclip_offsets(bc_num, platform, is_na_bc)
        row[f"Soft-clip Left clip pos offset vs Real BP Start ({plat_label})"]  = left_off
        row[f"Soft-clip Right clip pos offset vs Real BP End ({plat_label})"]   = right_off

    rows.append(row)

# ── Write output ──────────────────────────────────────────────────────────────
master = pd.DataFrame(rows)
master.to_csv(OUTPUT_TSV, sep="\t", index=False)
print(f"\nWritten: {OUTPUT_TSV}")
print(f"Shape  : {master.shape[0]} rows × {master.shape[1]} columns")
