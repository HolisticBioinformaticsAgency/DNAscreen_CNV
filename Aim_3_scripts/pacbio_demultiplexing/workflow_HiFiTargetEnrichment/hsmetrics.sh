#!/bin/bash
#SBATCH --job-name=hsmetrics_minimap2
#SBATCH --output=hsmetrics_minimap2_%j.out
#SBATCH --error=hsmetrics_minimap2_%j.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

### made this script because i don't want to have to run the snakemake workflow and mess up all my files ###

set -euo pipefail

############################
# Load modules
############################
module load picard
module load samtools

############################
# Paths
############################
BASE_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/demultiplex_pb/workflow_HiFiTargetEnrichment/batches/target_hifit"
REF="/fs04/vh83/reference/genomes/hg38/heng_li_recomended/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
REF_DICT="${REF%.fna}.dict"

# !! Set this to your enrichment panel BED file !!
TARGET_BED="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/3539131_Covered_DNA_Screen_for_ONT_alignment.bed"
INTERVAL_LIST="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/3539131_Covered_DNA_Screen_for_ONT_alignment.interval_list"

MULTIQC_OUT="${BASE_DIR}/multiqc_report"

############################
# Barcode list
############################
BARCODES=(
  bc1001 bc1002 bc1003 bc1004 bc1005 bc1006 bc1007 bc1008
  bc1097 bc1098 bc1099 bc1100 bc1101 bc1102 bc1103 bc1104
)

############################
# Step 1: Create sequence dictionary (once)
############################
if [[ ! -f "${REF_DICT}" ]]; then
  echo "Creating sequence dictionary..."
  picard CreateSequenceDictionary \
    R="${REF}" \
    O="${REF_DICT}"
fi

############################
# Step 2: BED -> interval list (once)
############################
if [[ ! -f "${INTERVAL_LIST}" ]]; then
  echo "Converting BED to interval list..."
  picard BedToIntervalList \
    I="${TARGET_BED}" \
    O="${INTERVAL_LIST}" \
    SD="${REF_DICT}" \
    SORT=true
fi

############################
# Step 3: CollectHsMetrics per barcode
############################
for BC in "${BARCODES[@]}"; do
  echo "===================================="
  echo "HsMetrics: ${BC}"
  echo "===================================="

  REALIGNED_BAM="${BASE_DIR}/${BC}/realigned_minimap2/${BC}.minimap2.GRCh38.bam"
  STATS_DIR="${BASE_DIR}/${BC}/realigned_minimap2_stats"

  mkdir -p "${STATS_DIR}"

  if [[ ! -f "${REALIGNED_BAM}" ]]; then
    echo "[WARNING] BAM not found for ${BC}, skipping"
    continue
  fi

  picard CollectHsMetrics \
    I="${REALIGNED_BAM}" \
    O="${STATS_DIR}/${BC}.hsmetrics.txt" \
    R="${REF}" \
    BAIT_INTERVALS="${INTERVAL_LIST}" \
    TARGET_INTERVALS="${INTERVAL_LIST}" \
    COVERAGE_CAP=1000 \
    NEAR_DISTANCE=250 \
    VALIDATION_STRINGENCY=LENIENT \
    TMP_DIR="${STATS_DIR}/tmp"

  echo "[${BC}] HsMetrics done"
done

############################
# Step 4: MultiQC (picks up .hsmetrics.txt automatically)
############################
echo "Running MultiQC..."
module load multiqc

mkdir -p "${MULTIQC_OUT}"

multiqc \
  --outdir "${MULTIQC_OUT}" \
  --filename "minimap2_multiqc_report" \
  --title "minimap2 HiFi Realignment QC" \
  "${BASE_DIR}"

echo "Done 🎉 Report: ${MULTIQC_OUT}/minimap2_multiqc_report.html"
