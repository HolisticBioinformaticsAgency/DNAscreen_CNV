#!/usr/bin/env python3



import os
import pandas as pd



# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = "/path/to/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers"
OUTPUT_TSV = os.path.join(BASE, "master_cnv_comparison_v3.tsv")



SOFTCLIP_RESULTS_DIR = (
    "/path/to/"
    "ONT_PacBio_CNV_calling/automated_precise_breakpoint_finding/"
    "batch_scripts/results/sv_analysis_with_sc_adjustment"
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
    "bc1113", "bc1114", "bc1115", "bc1116",
    "bc1117", "bc1118", "bc1119", "bc1125",
    "bc1126"
]




# ── aCGH reference coordinates (GRCh38, chrom without 'chr') ─────────────────
ACGH_COORDS = {
    "DNSxxxxxx": ("13", 32325998, 32326724),
    "DNSxxxxxx": ("3",  37046499, 37050922),
    "DNSxxxxxx": ("17", 43111887, 43117024),
    "DNSxxxxxx": ("2",  20900554, 21411082),
    "DNSxxxxxx": ("17", 43045165, 43045883),
    "DNSxxxxxx": ("16", 23606156, 23614296),
    "DNSxxxxxx": ("1",  55037565, 55040505),
    "DNSxxxxxx": ("17", 43078723, 43091151),
    "DNSxxxxxx": ("19", 11096325, 11166174),
    "DNSxxxxxx": ("2",  20997587, 21027015),
    # Inconclusive aCGH — coordinates exist but result is flagged as inconclusive
    "DNSxxxxxx": ("1",  55040596, 55060215),
    "DNSxxxxxx": ("3",  36993256, 37012096),
}



# Samples whose concordance column should always read "Inconclusive aCGH data"
INCONCLUSIVE_ACGH = {"DNSxxxxxx", "DNSxxxxxx"}



# For these samples, recompute offsets against the aCGH coords rather than the
# gene-window true_start/true_end stored in the comparison TSV. This takes
# priority over the generic gene-window fallback below, since aCGH coords are
# a tighter reference than the whole-gene body for these two samples.
OFFSET_OVERRIDE = {
    "DNSxxxxxx": ("1", 55040596, 55060215),
    "DNSxxxxxx": ("3", 36993256, 37012096),
}




# ── Load real breakpoints TSV ─────────────────────────────────────────────────
# Columns: Barcode | Sample ID | CNV event | Gene | True CNV Start |
#          True CNV End | Gene-window Coordinates | Notes
# True CNV Start/End format: "chrom:pos,pos,pos" (commas as thousand separators)
# Gene-window Coordinates format: "chrom:start,start-end,end" (commas as
# thousand separators within each number, hyphen separating start/end)
# This file is the SOURCE OF TRUTH for barcode -> Sample ID (DNS ID) mapping
# AND for whether a barcode has a known CNV event at all.
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




def parse_gene_window(coord_str):
    """Parse 'chrom:start-end' where start/end may contain comma
    thousand-separators. Returns (chrom_str, int_start, int_end) or
    (None, None, None) on failure / NA.
    """
    s = str(coord_str).strip()
    if s in ("NA", "", "nan"):
        return None, None, None
    try:
        chrom, rest = s.split(":", 1)
        start_str, end_str = rest.split("-", 1)
        start = int(start_str.replace(",", ""))
        end   = int(end_str.replace(",", ""))
        return chrom.strip(), start, end
    except (ValueError, AttributeError):
        return None, None, None




real_bp = {}  # barcode_num (str) -> {"sample_id": ..., "cnv_event": ..., "gene": ..., "true_start_*", "true_end_*", "gene_window_*"}
if os.path.exists(REAL_BREAKPOINTS_TSV):
    bp_df = pd.read_csv(REAL_BREAKPOINTS_TSV, sep="\t", dtype=str)
    bp_df.columns = bp_df.columns.str.strip()
    for _, bprow in bp_df.iterrows():
        bc_num = str(bprow.get("Barcode", "")).strip()
        sample_id = str(bprow.get("Sample ID", "")).strip()
        cnv_event = str(bprow.get("CNV event", "")).strip()
        gene = str(bprow.get("Gene", "")).strip()
        ts_chrom, ts_pos = parse_bp_coord(bprow.get("True CNV Start", "NA"))
        te_chrom, te_pos = parse_bp_coord(bprow.get("True CNV End", "NA"))
        gw_chrom, gw_start, gw_end = parse_gene_window(bprow.get("Gene-window Coordinates", "NA"))
        gw_raw = str(bprow.get("Gene-window Coordinates", "NA")).strip()
        real_bp[bc_num] = {
            "sample_id":         sample_id if sample_id not in ("", "nan") else "NA",
            "cnv_event":         cnv_event if cnv_event not in ("", "nan") else "NA",
            "gene":              gene if gene not in ("", "nan") else "NA",
            "true_start_chrom":  ts_chrom,
            "true_start_pos":    ts_pos,
            "true_end_chrom":    te_chrom,
            "true_end_pos":      te_pos,
            "gene_window_chrom": gw_chrom,
            "gene_window_start": gw_start,
            "gene_window_end":   gw_end,
            "gene_window_raw":   gw_raw if gw_raw not in ("", "nan") else "NA",
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




def compute_offsets(row, sample_id, note, is_na_barcode, bc_num):
    """
    Return (start_offset, end_offset) as formatted values.


    Priority of reference coordinates:
      1. OFFSET_OVERRIDE (aCGH coords) — for samples where aCGH is a tighter
         reference than the gene body (DNSxxxxxx, DNSxxxxxx).
      2. True breakpoints from real_breakpoints.tsv, via the tool TSV's
         precomputed start_offset/end_offset columns.
      3. Gene-window coordinates from real_breakpoints.tsv, when no true
         breakpoints exist for this barcode — offsets are recomputed here
         directly against the gene-window bounds and tagged accordingly.


    is_na_barcode reflects whether real_breakpoints.tsv has ANY known CNV
    event for this barcode (source of truth) — not the per-tool TSV's own
    note field, which may be stale relative to updates made directly to
    real_breakpoints.tsv.
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



    bp_info = real_bp.get(bc_num, {})
    has_true_bp = (
        bp_info.get("true_start_pos") is not None
        and bp_info.get("true_end_pos") is not None
    )



    if not has_true_bp:
        gw_start = bp_info.get("gene_window_start")
        gw_end   = bp_info.get("gene_window_end")
        if gw_start is not None and gw_end is not None:
            try:
                s_off = int(float(start)) - gw_start
                e_off = int(float(end))   - gw_end
                return f"{s_off} (vs gene-window)", f"{e_off} (vs gene-window)"
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




def get_concordant(row, sample_id, is_na_barcode):
    """
    is_na_barcode is resolved from real_breakpoints.tsv (source of truth),
    not from the per-tool TSV's own note field — see compute_offsets docstring.
    """
    note  = str(row.get("note", ""))
    start = row.get("best_match_start", "NA")



    if is_na_barcode:
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
    against a reference position for this barcode.


    Reference priority:
      1. True breakpoints (true_start_pos / true_end_pos) when available.
      2. Gene-window coordinates when no true breakpoints exist — offsets
         are tagged "(vs gene-window)" to flag the coarser reference.


    Returns ("N/A", "N/A") for NA barcodes or when neither reference exists.
    """
    if is_na_barcode:
        return "N/A", "N/A"



    bp_info = real_bp.get(bc_num)
    if bp_info is None:
        return "No real BP data", "No real BP data"



    has_true_bp = (
        bp_info["true_start_pos"] is not None
        and bp_info["true_end_pos"] is not None
    )
    has_gene_window = (
        bp_info.get("gene_window_start") is not None
        and bp_info.get("gene_window_end") is not None
    )



    if not has_true_bp and not has_gene_window:
        return "N/A", "N/A"



    if has_true_bp:
        ref_start, ref_end, ref_tag = bp_info["true_start_pos"], bp_info["true_end_pos"], ""
    else:
        ref_start, ref_end, ref_tag = bp_info["gene_window_start"], bp_info["gene_window_end"], " (vs gene-window)"



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
            off = int(float(str(clip_val).replace(",", ""))) - ref_pos
            return f"{off}{ref_tag}"
        except (ValueError, TypeError):
            return "Parse error"



    left_off  = clip_offset(left_raw,  ref_start)
    right_off = clip_offset(right_raw, ref_end)
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




# ── Get base info: real_breakpoints.tsv is the source of truth ───────────────
def get_base_info(bc):
    bc_num = bc.replace("bc", "")
    bp_info = real_bp.get(bc_num)



    if bp_info is not None:
        return {
            "Barcode":                 bc,
            "Sample ID":               bp_info["sample_id"],
            "CNV Type":                bp_info["cnv_event"],
            "Gene":                    bp_info.get("gene", "NA"),
            "True Start":              bp_info["true_start_pos"] if bp_info["true_start_pos"] is not None else "NA",
            "True End":                bp_info["true_end_pos"]   if bp_info["true_end_pos"]   is not None else "NA",
            "Gene-window Coordinates": bp_info.get("gene_window_raw", "NA"),
        }



    # Fallback (should rarely trigger): pull from whichever tool TSV has this
    # barcode, in case a barcode is missing from real_breakpoints.tsv entirely.
    for tool in TOOLS:
        for platform in PLATFORMS:
            df = data.get((tool, platform))
            if df is not None and bc in df.index:
                r = df.loc[bc].to_dict()
                return {
                    "Barcode":                 bc,
                    "Sample ID":               r.get("sample_id",  "NA"),
                    "CNV Type":                r.get("cnv_type",   "NA"),
                    "Gene":                    "NA",
                    "True Start":              r.get("true_start",  "NA"),
                    "True End":                r.get("true_end",    "NA"),
                    "Gene-window Coordinates": "NA",
                }
    return {"Barcode": bc, "Sample ID": "NA", "CNV Type": "NA", "Gene": "NA",
            "True Start": "NA", "True End": "NA", "Gene-window Coordinates": "NA"}




# ── Build master table ────────────────────────────────────────────────────────
rows = []



for bc in BARCODE_ORDER:
    base      = get_base_info(bc)
    row       = dict(base)
    sample_id = base["Sample ID"]
    bc_num    = bc.replace("bc", "")   # numeric portion, e.g. "1001"



    # Determine "has a known CNV event" status from real_breakpoints.tsv — the
    # source of truth — rather than from each per-tool TSV's own "note" column.
    # A per-tool TSV's note can be stale (generated before real_breakpoints.tsv
    # was updated with gene-window-only calls), which previously caused these
    # barcodes to be forced to "N/A" for coordinates even when a tool did call
    # something. Now a barcode is only treated as "no true CNV" when
    # real_breakpoints.tsv itself has no CNV event recorded for it.
    bc_cnv_event = real_bp.get(bc_num, {}).get("cnv_event", "NA")
    is_na_bc     = bc_cnv_event in ("NA", "nan")



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



            r    = df.loc[bc].to_dict()
            note = str(r.get("note", ""))



            s_off, e_off = compute_offsets(r, sample_id, note, is_na_bc, bc_num)



            row[col_coords]  = "N/A" if is_na_bc else format_coords(r)
            row[col_start]   = s_off
            row[col_end]     = e_off
            row[col_concord] = get_concordant(r, sample_id, is_na_bc)



    # ── Soft-clip offset columns (one pair per platform) ──────────────────────
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