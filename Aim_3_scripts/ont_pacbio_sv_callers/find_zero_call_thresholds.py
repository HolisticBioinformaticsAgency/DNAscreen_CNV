#!/usr/bin/env python3
"""
find_zero_call_thresholds_shared_platforms.py


For the copy-neutral male controls (bc1001, bc1099), sweep min SV size and
min read support thresholds to find the minimum combination that eliminates
all "doubly shared" SV calls overlapping the 9-gene BED, per tool.


A call is "doubly shared" if:
  1. It is SHARED between bc1001 and bc1099 within the SAME platform (ONT or
     PacBio) — same chromosome, start positions within POS_TOL bp, reciprocal
     size overlap ratio >= SIZE_RATIO_MIN.
  2. That barcode-shared call is ALSO shared between the ONT barcode-shared
     set and the PacBio barcode-shared set for that tool, using the same
     matching criteria (same chromosome, start within POS_TOL bp, reciprocal
     size ratio >= SIZE_RATIO_MIN).


In other words: first collapse bc1001+bc1099 -> barcode-shared calls (one set
per tool/platform), then collapse ONT+PacBio barcode-shared calls -> calls
shared across BOTH platforms (one set per tool). This isolates calls that are
reproducible across BOTH controls AND BOTH sequencing platforms — the most
likely systematic artifacts / germline variants.


Matching is done greedily, one-to-one, at each collapsing step. Size and
read_support are always taken as the max of the two matched values being
collapsed at that step.
"""


import os
import gzip
import re
import pandas as pd



# ── Paths ─────────────────────────────────────────────────────────────────────
BASE     = "/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers"
BED_PATH = os.path.join(BASE, "3539131_Covered_DNA_Screen_9genes.bed")


TOOLS     = ["sniffles", "cutesv", "svim", "debreak"]
PLATFORMS = ["ont", "pacbio"]
BARCODES  = ["bc1001", "bc1099"]


TOOL_LABELS = {
    "sniffles": "Sniffles",
    "cutesv":   "CuteSV",
    "svim":     "SVIM",
    "debreak":  "DeBreak",
}


# Matching tolerances for deciding whether two calls are the "same"
# underlying SV / artifact (used both for barcode-level and platform-level
# collapsing steps)
POS_TOL        = 100   # bp tolerance on start position
SIZE_RATIO_MIN = 0.7   # reciprocal size overlap ratio threshold (min(len)/max(len))


# Sweep ranges — adjust upper bounds if needed
SIZE_THRESHOLDS    = list(range(0, 50001, 50))   # 0, 50, 100, ..., 50000 bp
SUPPORT_THRESHOLDS = list(range(0, 101, 5))        # 0, 5, 10, ..., 100 reads



# ── Helpers (identical to main script) ───────────────────────────────────────
def normalise_chrom(c):
    return str(c).strip().lstrip("chr") if c else c



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



def parse_bed_overlapping_svs(vcf_path, tool, bed_regions):
    """Return list of dicts: chrom, start, end, size, support for BED-overlapping SVs."""
    svs = []
    if vcf_path is None:
        return svs
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


            rs = parse_read_support(info, fmt, samp, tool)
            svs.append({
                "chrom":   normalise_chrom(chrom),
                "start":   pos,
                "end":     sv_end,
                "size":    sv_end - pos,
                "support": rs if rs is not None else 0,
            })
    return svs



def find_shared_pairs(svs1, svs2, pos_tol=POS_TOL, size_ratio_min=SIZE_RATIO_MIN):
    """Greedy one-to-one match: return list of (sv1, sv2) pairs where sv1 is
    from svs1 and sv2 is from svs2, matched on same chrom, within pos_tol bp
    of start, with reciprocal size ratio >= size_ratio_min."""
    used2 = set()
    pairs = []
    for sv1 in svs1:
        best_j, best_dist = None, None
        for j, sv2 in enumerate(svs2):
            if j in used2 or sv1["chrom"] != sv2["chrom"]:
                continue
            dist = abs(sv1["start"] - sv2["start"])
            if dist > pos_tol:
                continue
            s1, s2 = max(sv1["size"], 1), max(sv2["size"], 1)
            ratio = min(s1, s2) / max(s1, s2)
            if ratio < size_ratio_min:
                continue
            if best_dist is None or dist < best_dist:
                best_dist, best_j = dist, j
        if best_j is not None:
            used2.add(best_j)
            pairs.append((sv1, svs2[best_j]))
    return pairs



def collapse_pairs(pairs):
    """Collapse matched (sv1, sv2) pairs into single records, keeping
    representative chrom/start (from sv1) and size/support = max of the two."""
    collapsed = []
    for sv1, sv2 in pairs:
        collapsed.append({
            "chrom":   sv1["chrom"],
            "start":   sv1["start"],
            "size":    max(sv1["size"], sv2["size"]),
            "support": max(sv1["support"], sv2["support"]),
        })
    return collapsed



# ── Step 1: Parse all BED-overlapping SVs per (barcode, tool, platform) ─────
bed_regions = load_bed(BED_PATH)
print(f"Loaded {len(bed_regions)} BED regions.\n")


raw_svs = {}
for barcode in BARCODES:
    for tool in TOOLS:
        for platform in PLATFORMS:
            key      = (barcode, tool, platform)
            vcf_path = find_vcf(tool, platform, barcode)
            if vcf_path is None:
                print(f"  WARNING: VCF not found — {barcode} | {tool} | {platform}")
                raw_svs[key] = []
                continue
            try:
                raw_svs[key] = parse_bed_overlapping_svs(vcf_path, tool, bed_regions)
            except Exception as e:
                print(f"  ERROR reading {vcf_path}: {e}")
                raw_svs[key] = []
            print(f"  Parsed {len(raw_svs[key]):4d} BED-overlapping SVs — "
                  f"{barcode} | {tool} | {platform}")



# ── Step 2: Collapse bc1001 vs bc1099 -> barcode-shared calls, per (tool, platform)
print("\nCollapsing bc1001 vs bc1099 (barcode-shared calls) per platform...\n")


barcode_shared = {}   # key: (tool, platform) -> list of collapsed dicts
for tool in TOOLS:
    for platform in PLATFORMS:
        svs_1001 = raw_svs.get(("bc1001", tool, platform), [])
        svs_1099 = raw_svs.get(("bc1099", tool, platform), [])
        pairs    = find_shared_pairs(svs_1001, svs_1099)
        barcode_shared[(tool, platform)] = collapse_pairs(pairs)
        print(f"  {TOOL_LABELS[tool]:8s} | {platform.upper():6s} → "
              f"{len(pairs)} barcode-shared calls "
              f"(bc1001={len(svs_1001)}, bc1099={len(svs_1099)})")



# ── Step 3: Collapse ONT vs PacBio barcode-shared calls -> doubly-shared, per tool
print("\nCollapsing ONT vs PacBio barcode-shared calls (doubly-shared) per tool...\n")


doubly_shared = {}   # key: tool -> list of {size, support} dicts
for tool in TOOLS:
    ont_shared    = barcode_shared.get((tool, "ont"), [])
    pacbio_shared = barcode_shared.get((tool, "pacbio"), [])
    pairs         = find_shared_pairs(ont_shared, pacbio_shared)
    doubly_shared[tool] = [
        {"size": max(sv1["size"], sv2["size"]),
         "support": max(sv1["support"], sv2["support"])}
        for sv1, sv2 in pairs
    ]
    print(f"  {TOOL_LABELS[tool]:8s} → {len(pairs)} doubly-shared calls "
          f"(ONT barcode-shared={len(ont_shared)}, "
          f"PacBio barcode-shared={len(pacbio_shared)})")



# ── Sweep thresholds on doubly-shared calls only ─────────────────────────────
print("\nSweeping thresholds on doubly-shared calls (barcode + platform)...\n")


results = []


for min_size in SIZE_THRESHOLDS:
    for min_support in SUPPORT_THRESHOLDS:
        total_remaining = 0
        breakdown = {}
        for tool in TOOLS:
            svs  = doubly_shared.get(tool, [])
            kept = [sv for sv in svs
                    if sv["size"] >= min_size
                    and sv["support"] >= min_support]
            breakdown[tool] = len(kept)
            total_remaining += len(kept)


        results.append({
            "min_size_bp":      min_size,
            "min_support":      min_support,
            "total_remaining":  total_remaining,
            "breakdown":        breakdown,
        })


        if total_remaining == 0:
            print(f"  ZERO doubly-shared calls reached → min_size >= {min_size} bp, "
                  f"min_support >= {min_support} reads")
            break
    else:
        continue
    break



# ── Report ────────────────────────────────────────────────────────────────────
zero_results = [r for r in results if r["total_remaining"] == 0]


if not zero_results:
    print("\nNo threshold combination within the sweep range achieved zero doubly-shared calls.")
    print("Consider increasing SIZE_THRESHOLDS or SUPPORT_THRESHOLDS upper bounds.\n")
else:
    best = zero_results[0]
    print(f"\n{'='*60}")
    print(f"Minimum thresholds for zero doubly-shared BED-overlapping calls:")
    print(f"  Min SV size    : >= {best['min_size_bp']} bp")
    print(f"  Min read support: >= {best['min_support']} reads")
    print(f"{'='*60}\n")
    print("Per-tool doubly-shared call counts at these thresholds:")
    for tool, count in sorted(best["breakdown"].items()):
        print(f"  {TOOL_LABELS[tool]:8s} → {count} doubly-shared calls")


# ── Write full sweep results to TSV ──────────────────────────────────────────
OUTPUT_TSV = os.path.join(BASE, "threshold_sweep_results_shared_platforms.tsv")
sweep_rows = []
for r in results:
    row = {"min_size_bp": r["min_size_bp"], "min_support": r["min_support"],
           "total_remaining": r["total_remaining"]}
    for tool, count in r["breakdown"].items():
        row[f"{tool}"] = count
    sweep_rows.append(row)


sweep_df = pd.DataFrame(sweep_rows)
sweep_df.to_csv(OUTPUT_TSV, sep="\t", index=False)
print(f"\nFull sweep written to: {OUTPUT_TSV}")
