#!/bin/bash
#SBATCH --job-name=sa_junction_counts
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/path/to/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/sa_junction_counts_%j.out
#SBATCH --error=/path/to/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/sa_junction_counts_%j.err

set -euo pipefail

source /usr/local/anaconda/5.1.0-Python3.6-gcc5/etc/profile.d/conda.sh
conda activate ont
module load samtools

BASE="/path/to/ONT_PacBio_CNV_calling"
BASE_ONT="${BASE}/breakpoint_identification_ont"
BASE_PB="${BASE}/breakpoint_identification_pacbio"

# ── ONT full-BAM source directories ────────────────────────────────────────
# Barcodes 1113-1128 come from a different dorado demultiplex run than the
# rest of the cohort (1001-1104), so pick the base dir per barcode.
FULL_BAM_ONT_DIR_DEFAULT="/path/to/dorado_demultiplex_142pods_barcode_both_ends/bams_dedup_cleaned"
FULL_BAM_ONT_DIR_ROUND2="/path/to/dorado_demultiplex_round2_full_ont"

FULL_BAM_PB_DIR="/path/to/demultiplex_pb/workflow_HiFiTargetEnrichment/batches/target_hifit"
OUTPUT="${BASE}/ont_pacbio_sv_callers/sa_junction_read_counts_v3.tsv"

get_full_bam_ont_dir() {
    local bc="$1"
    if (( bc >= 1113 && bc <= 1128 )); then
        echo "$FULL_BAM_ONT_DIR_ROUND2"
    else
        echo "$FULL_BAM_ONT_DIR_DEFAULT"
    fi
}

# ── N/A barcodes ──────────────────────────────────────────────────────────────
declare -A NA_SAMPLE_ID
NA_SAMPLE_ID[1001]="Copy Neutral Male"
NA_SAMPLE_ID[1006]="DNSxxxxxx"
NA_SAMPLE_ID[1099]="Copy Neutral Male"
NA_SAMPLE_ID[1103]="DNSxxxxxx"

# ── Junction sequences ────────────────────────────────────────────────────────
declare -A JUNCTION
declare -A SAMPLE_ID

JUNCTION[1002]="gctaatttttg";                                                    SAMPLE_ID[1002]="DNSxxxxxx"
JUNCTION[1003]="ctgcccaccttggcctcccaaagtgctgggattacaggcgtgagccactgcacctggcc";   SAMPLE_ID[1003]="DNSxxxxxx"
JUNCTION[1004]="caggagatcgagaccagcctggccaa";                                     SAMPLE_ID[1004]="DNSxxxxxx"
JUNCTION[1005]="cacgccattctcactaggatctag";                                        SAMPLE_ID[1005]="DNSxxxxxx"
JUNCTION[1007]="tgggaacagcagaggtttctttgtt";                                      SAMPLE_ID[1007]="DNSxxxxxx"
JUNCTION[1008]="tttcaccatgttggccaggctggtcttaaactc";                              SAMPLE_ID[1008]="DNSxxxxxx"
JUNCTION[1097]="tggccctttcacctgaggctccgagaggtt";                                 SAMPLE_ID[1097]="DNSxxxxxx"
JUNCTION[1098]="caggagatcgagaccagcctggccaa";                                     SAMPLE_ID[1098]="DNSxxxxxx"
JUNCTION[1100]="ctgcccaccttggcctcccaaagtgctgggattacaggcgtgagccactgcacctggcc";   SAMPLE_ID[1100]="DNSxxxxxx"
JUNCTION[1101]="tttcaccatgttggccaggctggtcttaaactc";                              SAMPLE_ID[1101]="DNSxxxxxx"
JUNCTION[1102]="tttcaccatgttggccaggctggtcttaaactc";                              SAMPLE_ID[1102]="DNSxxxxxx"
JUNCTION[1104]="aaatggggtctcctgcctca";                                            SAMPLE_ID[1104]="DNSxxxxxx"
JUNCTION[1113]="ctcggcctcccaaagtgctgggattacaggcgtga";                            SAMPLE_ID[1113]="DNSxxxxxx"
JUNCTION[1114]="ctctgtctgtttc";                                                   SAMPLE_ID[1114]="DNSxxxxxx"

# TODO: fill in JUNCTION / SAMPLE_ID for barcodes 1115-1128 here, e.g.:
# JUNCTION[1115]="..."; SAMPLE_ID[1115]="..."
# ...
# JUNCTION[1128]="..."; SAMPLE_ID[1128]="..."

# ── True breakpoint coordinates (chrom without commas) ───────────────────────
declare -A BP_CHROM BP_START BP_END

BP_CHROM[1002]="16"; BP_START[1002]="23605575";  BP_END[1002]="23615114"
BP_CHROM[1003]="17"; BP_START[1003]="43045048";  BP_END[1003]="43046166"
BP_CHROM[1004]="17"; BP_START[1004]="43111616";  BP_END[1004]="43117534"
BP_CHROM[1005]="13"; BP_START[1005]="32325651";  BP_END[1005]="32327136"
BP_CHROM[1007]="2";  BP_START[1007]="20993431";  BP_END[1007]="21026997"
BP_CHROM[1008]="17"; BP_START[1008]="43078282";  BP_END[1008]="43084362"
BP_CHROM[1097]="1";  BP_START[1097]="55029215";  BP_END[1097]="55043368"
BP_CHROM[1098]="17"; BP_START[1098]="43111616";  BP_END[1098]="43117534"
BP_CHROM[1100]="17"; BP_START[1100]="43045048";  BP_END[1100]="43046166"
BP_CHROM[1101]="17"; BP_START[1101]="43078282";  BP_END[1101]="43084362"
BP_CHROM[1102]="17"; BP_START[1102]="43078282";  BP_END[1102]="43084362"
BP_CHROM[1104]="19"; BP_START[1104]="11094251";  BP_END[1104]="11175478"
# PALB2 tandem duplication
BP_CHROM[1113]="16"; BP_START[1113]="23611399";  BP_END[1113]="23615201"
# MLH1 tandem duplication
BP_CHROM[1114]="3";  BP_START[1114]="37045688";  BP_END[1114]="37052535"

# TODO: fill in BP_CHROM / BP_START / BP_END for barcodes 1115-1128 here.

# ── Helpers ───────────────────────────────────────────────────────────────────
revcomp() {
    echo "$1" | tr 'acgtACGT' 'tgcaTGCA' | rev
}

count_junction_reads() {
    local bam="$1" fwd="$2" rc="$3"
    samtools view "$bam" \
    | awk -v fwd="$fwd" -v rc="$rc" '{
        seq = tolower($10)
        if (index(seq, fwd) > 0 || index(seq, rc) > 0) print $1
    }' | sort -u | wc -l
}

count_locus_reads() {
    local bam="$1" chrom="$2" bp1="$3" bp2="$4"
    {
        samtools view "$bam" "chr${chrom}:${bp1}-${bp1}"
        samtools view "$bam" "chr${chrom}:${bp2}-${bp2}"
    } | awk '{print $1}' | sort -u | wc -l
}

pct() {
    local num="$1" denom="$2"
    if [[ "$denom" -eq 0 ]]; then echo "N/A"; return; fi
    awk -v n="$num" -v d="$denom" 'BEGIN { printf "%.1f%%", n*100/d }'
}

find_bam() {
    local pattern="$1"
    local hit
    hit=$(ls ${pattern} 2>/dev/null | head -1 || true)
    echo "${hit:-NOT_FOUND}"
}

# ── Header ────────────────────────────────────────────────────────────────────
echo -e "barcode\tsample_id\tjunction_seq\t\
SA_junction_reads_ONT\ttotal_reads_at_breakpoints_ONT\tpct_junction_reads_ONT\t\
SA_junction_reads_PacBio\ttotal_reads_at_breakpoints_PacBio\tpct_junction_reads_PacBio\t\
full_bam_ont\tbreakpoint_reads_bam_ont\t\
full_bam_pacbio\tbreakpoint_reads_bam_pacbio" > "$OUTPUT"

# ── Main loop ─────────────────────────────────────────────────────────────────
# Extended to 1128 to cover the round-2 ONT barcodes.
for bc_num in 1001 1002 1003 1004 1005 1006 1007 1008 \
              1097 1098 1099 1100 1101 1102 1103 1104 \
              1113 1114 1115 1116 1117 1118 1119 1120 \
              1121 1122 1123 1124 1125 1126 1127 1128; do

    # N/A barcodes ─────────────────────────────────────────────────────────────
    if [[ -n "${NA_SAMPLE_ID[$bc_num]+_}" ]]; then
        sid="${NA_SAMPLE_ID[$bc_num]}"
        echo "── barcode${bc_num} (${sid}) — N/A"
        echo -e "barcode${bc_num}\t${sid}\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A" >> "$OUTPUT"
        continue
    fi

    # Skip barcodes that have no junction/breakpoint data yet (e.g. undefined
    # 1115-1128 until you populate JUNCTION/SAMPLE_ID/BP_* above).
    if [[ -z "${JUNCTION[$bc_num]+_}" ]]; then
        echo "── barcode${bc_num} — SKIPPED (no junction/sample/breakpoint data defined)"
        echo -e "barcode${bc_num}\tUNDEFINED\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A" >> "$OUTPUT"
        continue
    fi

    junction="${JUNCTION[$bc_num]}"
    junction_rc=$(revcomp "$junction")
    sid="${SAMPLE_ID[$bc_num]}"
    chrom="${BP_CHROM[$bc_num]}"
    bp1="${BP_START[$bc_num]}"
    bp2="${BP_END[$bc_num]}"

    echo "── barcode${bc_num} (${sid})  breakpoints: chr${chrom}:${bp1} / chr${chrom}:${bp2}"

    # Full BAM paths ───────────────────────────────────────────────────────────
    full_bam_ont_dir=$(get_full_bam_ont_dir "$bc_num")
    full_bam_ont="${full_bam_ont_dir}/barcode${bc_num}_q10_mq20_dedup.bam"
    [[ ! -f "$full_bam_ont" ]] && full_bam_ont="NOT_FOUND"

    full_bam_pb="${FULL_BAM_PB_DIR}/bc${bc_num}/realigned_minimap2/bc${bc_num}.minimap2.GRCh38.bam"
    [[ ! -f "$full_bam_pb" ]] && full_bam_pb="NOT_FOUND"

    # Breakpoint + locus BAM paths ─────────────────────────────────────────────
    bp_bam_ont=$(find_bam   "${BASE_ONT}/barcode${bc_num}/barcode${bc_num}_*.breakpoint_reads.bam")
    locus_bam_ont=$(find_bam "${BASE_ONT}/barcode${bc_num}/barcode${bc_num}_*.locus.bam")

    bp_bam_pb=$(find_bam    "${BASE_PB}/barcode${bc_num}/minimap2-aligned/bc${bc_num}.*.breakpoint_reads.bam")
    locus_bam_pb=$(find_bam  "${BASE_PB}/barcode${bc_num}/minimap2-aligned/bc${bc_num}.*.locus.bam")

    # ONT counts ───────────────────────────────────────────────────────────────
    if [[ "$bp_bam_ont" != "NOT_FOUND" ]]; then
        junc_ont=$(count_junction_reads "$bp_bam_ont" "$junction" "$junction_rc")
    else
        junc_ont="NO_BAM"
    fi

    if [[ "$locus_bam_ont" != "NOT_FOUND" ]]; then
        total_ont=$(count_locus_reads "$locus_bam_ont" "$chrom" "$bp1" "$bp2")
        if [[ "$junc_ont" == "NO_BAM" ]]; then
            pct_ont="NO_BAM"
        else
            pct_ont=$(pct "$junc_ont" "$total_ont")
        fi
    else
        total_ont="NO_BAM"; pct_ont="NO_BAM"
    fi

    echo "  [ONT] junction reads: ${junc_ont} / ${total_ont}  (${pct_ont})"

    # PacBio counts ────────────────────────────────────────────────────────────
    if [[ "$bp_bam_pb" != "NOT_FOUND" ]]; then
        junc_pb=$(count_junction_reads "$bp_bam_pb" "$junction" "$junction_rc")
    else
        junc_pb="NO_BAM"
    fi

    if [[ "$locus_bam_pb" != "NOT_FOUND" ]]; then
        total_pb=$(count_locus_reads "$locus_bam_pb" "$chrom" "$bp1" "$bp2")
        if [[ "$junc_pb" == "NO_BAM" ]]; then
            pct_pb="NO_BAM"
        else
            pct_pb=$(pct "$junc_pb" "$total_pb")
        fi
    else
        total_pb="NO_BAM"; pct_pb="NO_BAM"
    fi

    echo "  [PB]  junction reads: ${junc_pb} / ${total_pb}  (${pct_pb})"

    # Write row ────────────────────────────────────────────────────────────────
    echo -e "barcode${bc_num}\t${sid}\t${junction}\t\
${junc_ont}\t${total_ont}\t${pct_ont}\t\
${junc_pb}\t${total_pb}\t${pct_pb}\t\
${full_bam_ont}\t${bp_bam_ont}\t\
${full_bam_pb}\t${bp_bam_pb}" >> "$OUTPUT"

done

echo ""
echo "Done. Results written to: ${OUTPUT}"
