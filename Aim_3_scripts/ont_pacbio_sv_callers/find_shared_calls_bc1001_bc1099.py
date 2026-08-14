#!/usr/bin/env python3
"""
find_shared_calls_bc1001_bc1099.py


For the copy-neutral male controls (bc1001, bc1099), for each SV caller and
each platform (ONT / PacBio) separately, count how many BED-overlapping SV
calls are SHARED between bc1001 and bc1099 (i.e. likely germline/artifactual
calls common to both controls, rather than being unique to one control).


A call in bc1001 is considered "shared" with a call in bc1099 if:
  - same chromosome
  - start positions within POS_TOL bp of each other
  - reciprocal size overlap ratio >= SIZE_RATIO_MIN


Matching is done greedily, one-to-one (each bc1099 call can only be consumed
once), to avoid double counting when there are clusters of nearby calls.


For each shared call pair, the max read_support and max SV size (taken as
the larger of the two matched bc1001/bc1099 values) is recorded. The single
shared pair with the overall maximum read_support is reported (along with
that pair's SV size), and the single shared pair with the overall maximum
SV size is reported (along with that pair's read_support), per (tool, platform).


Output: shared_calls_bc1001_bc1099.tsv
Columns: tool, platform, n_bc1001_calls, n_bc1099_calls, n_shared_calls,
         frac_shared_of_bc1001, frac_shared_of_bc1099,
         max_read_support, max_read_support_sv_size,
         max_sv_size, max_sv_size_read_support
"""


import os
import gzip
import re
import pandas as pd


# ── Paths ─────────────────────────────────────────────────────────────────
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


# Matching tolerances for deciding whether two calls (one per barcode) are
# the "same" underlying SV / artifact
POS_TOL      = 100   # bp tolerance on start position
SIZE_RATIO_MIN = 0.7 # reciprocal size ratio threshold (min(len)/max(len))


OUTPUT_TSV = os.path.join(BASE, "shared_calls_bc1001_bc1099.tsv")



# ── Helpers (identical to main script) ──────────────────────────────────────
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
    from svs1 (bc1001) and sv2 is from svs2 (bc1099), matched on same chrom,
    within pos_tol bp of start, with reciprocal size ratio >= size_ratio_min."""
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



# ── Main ─────────────────────────────────────────────────────────────────
bed_regions = load_bed(BED_PATH)
print(f"Loaded {len(bed_regions)} BED regions.\n")


rows = []
for tool in TOOLS:
    for platform in PLATFORMS:
        svs_by_barcode = {}
        for barcode in BARCODES:
            vcf_path = find_vcf(tool, platform, barcode)
            if vcf_path is None:
                print(f"  WARNING: VCF not found — {barcode} | {tool} | {platform}")
                svs_by_barcode[barcode] = []
                continue
            try:
                svs_by_barcode[barcode] = parse_bed_overlapping_svs(vcf_path, tool, bed_regions)
            except Exception as e:
                print(f"  ERROR reading {vcf_path}: {e}")
                svs_by_barcode[barcode] = []


        svs_1001 = svs_by_barcode["bc1001"]
        svs_1099 = svs_by_barcode["bc1099"]
        pairs    = find_shared_pairs(svs_1001, svs_1099)
        n_shared = len(pairs)


        n1, n2 = len(svs_1001), len(svs_1099)
        frac1 = n_shared / n1 if n1 else 0.0
        frac2 = n_shared / n2 if n2 else 0.0


        # Find the shared pair with max read_support (report its SV size too)
        # and the shared pair with max SV size (report its read_support too).
        # Values taken as the larger of the two matched bc1001/bc1099 values.
        max_support_val, max_support_sv_size = None, None
        max_size_val, max_size_read_support  = None, None


        for sv1, sv2 in pairs:
            pair_max_support = max(sv1["support"], sv2["support"])
            pair_max_size    = max(sv1["size"], sv2["size"])


            if max_support_val is None or pair_max_support > max_support_val:
                max_support_val      = pair_max_support
                max_support_sv_size  = pair_max_size


            if max_size_val is None or pair_max_size > max_size_val:
                max_size_val          = pair_max_size
                max_size_read_support = pair_max_support


        rows.append({
            "tool":                       TOOL_LABELS[tool],
            "platform":                   platform.upper(),
            "n_bc1001_calls":             n1,
            "n_bc1099_calls":             n2,
            "n_shared_calls":             n_shared,
            "frac_shared_of_bc1001":      round(frac1, 4),
            "frac_shared_of_bc1099":      round(frac2, 4),
            "max_read_support":           max_support_val,
            "max_read_support_sv_size":   max_support_sv_size,
            "max_sv_size":                max_size_val,
            "max_sv_size_read_support":   max_size_read_support,
        })


        print(f"  {TOOL_LABELS[tool]:8s} | {platform.upper():6s} → "
              f"bc1001={n1}, bc1099={n2}, shared={n_shared}, "
              f"max_support={max_support_val} (size={max_support_sv_size}), "
              f"max_size={max_size_val} (support={max_size_read_support})")


shared_df = pd.DataFrame(rows)
shared_df.to_csv(OUTPUT_TSV, sep="\t", index=False)
print(f"\nShared-call summary written to: {OUTPUT_TSV}")