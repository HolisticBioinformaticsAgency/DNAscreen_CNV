#!/bin/bash
#SBATCH --job-name=dorado_sup_round2_full
#SBATCH --output=dorado_sup_round2_full_%j.out
#SBATCH --error=dorado_sup_round2_full_%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

module load samtools

# ---- Record total runtime ---- #
START_TIME=$(date +%s)
echo "🚀 Job started at $(date)"

format_duration() {
    local total_seconds=$1
    local hours=$((total_seconds / 3600))
    local minutes=$(((total_seconds % 3600) / 60))
    local seconds=$((total_seconds % 60))
    printf "%02dh:%02dm:%02ds" "${hours}" "${minutes}" "${seconds}"
}

report_runtime() {
    local exit_status=$?
    local end_time
    local total_seconds

    end_time=$(date +%s)
    total_seconds=$((end_time - START_TIME))

    echo "⏱️ Total processing time: $(format_duration "${total_seconds}")"
    echo "🏁 Job ended at $(date)"

    if [[ ${exit_status} -eq 0 ]]; then
        echo "✅ Pipeline completed successfully"
    else
        echo "❌ Pipeline failed with exit status ${exit_status}"
    fi

    exit "${exit_status}"
}

trap report_runtime EXIT

# ---- Define paths ---- #
DORADO_BIN="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_test/dorado-1.2.0-linux-x64/bin/dorado"
MODEL="/fs04/vh83/ont_dnascreen/jason_test/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
POD5_DIR="/fs04/vh83/ont_dnascreen/round2_pod5_last_40"
BASE_OUTPUT_DIR="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_round2_full_ont"
OUTPUT_DIR="${BASE_OUTPUT_DIR}/dorado_sup_output"
BAM_FILE="${OUTPUT_DIR}/PBG10946_pass_fb6074ec_f8186473_0_custom_barcode_0_19.bam"
ARRANGEMENT="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_round2_full_ont/test_arrangement_extended_mask.toml"
BARCODES="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_round2_full_ont/barcodes_RC_F.fa"
KIT_NAME="custom_barcode"
REFERENCE="/fs04/vh83/reference/genomes/hg38/heng_li_recomended/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# ---- Path to the first-87 round BAM (the base for concatenation) ---- #
FIRST87_BAM="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_round2_first87_ont/dorado_sup_output/PBG10946_pass_fb6074ec_f8186473_0_custom_barcode_0_19.bam"
COMBINED_BAM="${OUTPUT_DIR}/PBG10946_pass_fb6074ec_f8186473_0_custom_barcode_0_19_combined.bam"

# DEMUX_DIR is shared with rename_part.sh — keep this value identical
# in both scripts if you ever change it.
DEMUX_DIR="${BASE_OUTPUT_DIR}/bams_demuxed"

mkdir -p "${OUTPUT_DIR}"

# ---- Run Dorado ---- #
echo "⏳ Running Dorado on all POD5 files in ${POD5_DIR}"
DORADO_START=$(date +%s)

"${DORADO_BIN}" basecaller \
    "${MODEL}" \
    "${POD5_DIR}" \
    --barcode-arrangement "${ARRANGEMENT}" \
    --barcode-sequences "${BARCODES}" \
    --kit-name "${KIT_NAME}" \
    --reference "${REFERENCE}" \
    --barcode-both-ends \
    > "${BAM_FILE}"

DORADO_END=$(date +%s)
echo "✅ Dorado SUP basecalling completed at $(date)"
echo "⏱️ Dorado runtime: $(format_duration "$((DORADO_END - DORADO_START))")"
echo "📦 BAM saved as ${BAM_FILE}"

# ---- Sanity check on both BAMs before touching anything ---- #
if [[ ! -f "${BAM_FILE}" ]]; then
    echo "⚠️ Round-2 (last-40) BAM file not found: ${BAM_FILE}"
    exit 1
fi

if [[ ! -f "${FIRST87_BAM}" ]]; then
    echo "⚠️ First-87 BAM file not found: ${FIRST87_BAM}"
    exit 1
fi

# ---- Concatenate first-87 BAM (base) + last-40 BAM (on top) ---- #
echo "🔗 Concatenating BAMs: ${FIRST87_BAM} + ${BAM_FILE}"
CONCAT_START=$(date +%s)

# samtools cat requires the input BAMs to share a compatible header
# (same @SQ reference lines). It does NOT sort or re-index; it just
# appends records in the order given. -h takes the header from the
# first file listed.
samtools cat \
    -h "${FIRST87_BAM}" \
    -o "${COMBINED_BAM}" \
    "${FIRST87_BAM}" \
    "${BAM_FILE}"

CONCAT_END=$(date +%s)
echo "✅ Concatenation completed at $(date)"
echo "⏱️ Concat runtime: $(format_duration "$((CONCAT_END - CONCAT_START))")"
echo "📦 Combined BAM saved as ${COMBINED_BAM}"

# From this point on, splitting operates on the COMBINED BAM, not the
# round-2-only BAM_FILE, so downstream barcode bins contain reads from
# all 127 pod5 files (first 87 + last 40).
BAM_FILE="${COMBINED_BAM}"

# ---- Split BAM by barcode ---- #
if [[ -f "${BAM_FILE}" ]]; then
    SPLIT_START=$(date +%s)
    echo "📦 Splitting combined BAM by barcode and generating unaccounted.bam"

    mkdir -p "${DEMUX_DIR}"

    samtools split \
        -d BC:Z \
        "${BAM_FILE}" \
        -u "${DEMUX_DIR}/unaccounted.bam" \
        -f "${DEMUX_DIR}/PBG10946_pass_fb6074ec_f8186473_0_custom_barcode_%#.bam"

    SPLIT_END=$(date +%s)
    echo "⏱️ BAM split runtime: $(format_duration "$((SPLIT_END - SPLIT_START))")"
    echo "📂 Split BAMs written to ${DEMUX_DIR}"
    echo "➡️  Next step: submit rename_part.sh with --dependency=afterok:${SLURM_JOB_ID:-<this_jobid>}"
else
    echo "⚠️ BAM file not found: ${BAM_FILE}"
    exit 1
fi

echo "🎉 Dorado → concat → split stage completed successfully at $(date)"
