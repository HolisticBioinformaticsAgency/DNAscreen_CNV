#!/bin/bash
#SBATCH --job-name=dorado_sup_all
#SBATCH --output=dorado_sup_all_%j.out
#SBATCH --error=dorado_sup_all_%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

module load samtools

# ---- Define paths ---- #
DORADO_BIN="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_test/dorado-1.2.0-linux-x64/bin/dorado"
MODEL="/fs04/vh83/ont_dnascreen/jason_test/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
POD5_DIR="/fs04/vh83/ont_dnascreen/pod5"
BASE_OUTPUT_DIR="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends"
OUTPUT_DIR="${BASE_OUTPUT_DIR}/dorado_sup_output"
BAM_FILE="${OUTPUT_DIR}/PBG10946_pass_fb6074ec_f8186473_0_custom_barcode_0_19.bam"
ARRANGEMENT="/fs04/vh83/jason/ont/test_arrangement_extended_mask.toml"
BARCODES="/fs04/vh83/jason/ont/barcodes_RC_F.fa"
KIT_NAME="custom_barcode"
REFERENCE="/fs04/vh83/reference/genomes/hg38/heng_li_recomended/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

mkdir -p "$OUTPUT_DIR"

# # ---- Run Dorado ---- #
echo "⏳ Running Dorado on all pod5 files in ${POD5_DIR}"
"${DORADO_BIN}" basecaller \
    "${MODEL}" \
    "${POD5_DIR}" \
    --barcode-arrangement "${ARRANGEMENT}" \
    --barcode-sequences "${BARCODES}" \
    --kit-name "${KIT_NAME}" \
    --reference "${REFERENCE}" \
    --barcode-both-ends \
    > "${BAM_FILE}"

echo "✅ Dorado SUP basecalling completed at $(date)"
echo "📦 BAM saved as ${BAM_FILE}"

# ---- Split BAM by barcode, capture unaccounted reads ---- #
UNACCOUNTED_BAM="${OUTPUT_DIR}/unaccounted.bam"

if [[ -f "${BAM_FILE}" ]]; then
    echo "📦 Splitting BAM by barcode and generating unaccounted.bam"

    # Create output subdirectory for demuxed BAMs
    DEMUX_DIR="${OUTPUT_DIR}/bams_demuxed"
    mkdir -p "${DEMUX_DIR}"

    # Run samtools split
    samtools split -d BC:Z "${BAM_FILE}" -u "${DEMUX_DIR}/unaccounted.bam" -f "${DEMUX_DIR}/PBG10946_pass_fb6074ec_f8186473_0_custom_barcode_%#.bam"

    # Index all split BAMs
    for split_bam in "${DEMUX_DIR}"/PBG10946_pass_fb6074ec_f8186473_0_custom_barcode_*.bam "${DEMUX_DIR}/unaccounted.bam"; do
        if [[ -f "${split_bam}" ]]; then
            echo "📌 Indexing ${split_bam}"
            samtools index "${split_bam}"
        fi
    done
else
    echo "⚠️ BAM file not found: ${BAM_FILE}"
fi

echo "🎉 All processing completed successfully at $(date)"

