#!/bin/bash
#SBATCH --job-name=cutesv
#SBATCH --array=0-15
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/cutesv/logs/%A_%a.out
#SBATCH --error=/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/cutesv/logs/%A_%a.err

set -euo pipefail

# ── Set mode: "ont" or "pacbio" ──────────────────────────────────────────────
MODE="ont"   # <-- Change to "pacbio" for PacBio runs
             #     Also update --job-name above to "cutesv_pacbio"!
# ─────────────────────────────────────────────────────────────────────────────

source /usr/local/anaconda/5.1.0-Python3.6-gcc5/etc/profile.d/conda.sh
conda activate ont

REFERENCE="/fs04/vh83/reference/genomes/hg38/heng_li_recomended/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Full list of barcodes (indices 0-15)
BARCODES=(
    bc1001 bc1002 bc1003 bc1004 bc1005 bc1006 bc1007 bc1008
    bc1097 bc1098 bc1099 bc1100 bc1101 bc1102 bc1103 bc1104
)

barcode=${BARCODES[$SLURM_ARRAY_TASK_ID]}
THREADS=${SLURM_CPUS_PER_TASK}

if [[ "${MODE}" == "ont" ]]; then
    BAM_DIR="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/bams_dedup_cleaned"
    OUTBASE="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/cutesv/cutesv_ont"
    bc_num="${barcode#bc}"
    bam="${BAM_DIR}/barcode${bc_num}_q10_mq20_dedup.bam"
    MAX_BIAS_INS=100;  RATIO_INS=0.3
    MAX_BIAS_DEL=100;  RATIO_DEL=0.3

elif [[ "${MODE}" == "pacbio" ]]; then
    BASE="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/demultiplex_pb/workflow_HiFiTargetEnrichment/batches/target_hifit"
    OUTBASE="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/cutesv/cutesv_pacbio"
    bam="${BASE}/${barcode}/realigned_minimap2/${barcode}.minimap2.GRCh38.bam"
    MAX_BIAS_INS=1000; RATIO_INS=0.9
    MAX_BIAS_DEL=1000; RATIO_DEL=0.5

else
    echo "ERROR: MODE must be 'ont' or 'pacbio'. Got: ${MODE}" >&2
    exit 1
fi

vcf="${OUTBASE}/${barcode}/cutesv_${barcode}.vcf"
work_dir="${OUTBASE}/${barcode}/work"

mkdir -p "${OUTBASE}/logs"
mkdir -p "${OUTBASE}/${barcode}"
mkdir -p "${work_dir}"

if [[ ! -f "${bam}" ]]; then
    echo "ERROR: BAM not found: ${bam}" >&2
    exit 1
fi

echo "[$(date)] Running cuteSV (${MODE}) for ${barcode} (task ${SLURM_ARRAY_TASK_ID})"
cuteSV "${bam}" "${REFERENCE}" "${vcf}" "${work_dir}" \
    --threads "${THREADS}" \
    --max_cluster_bias_INS "${MAX_BIAS_INS}" \
    --diff_ratio_merging_INS "${RATIO_INS}" \
    --max_cluster_bias_DEL "${MAX_BIAS_DEL}" \
    --diff_ratio_merging_DEL "${RATIO_DEL}"
echo "[$(date)] Done: ${barcode}"