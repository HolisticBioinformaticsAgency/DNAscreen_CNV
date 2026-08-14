#!/bin/bash
#SBATCH --job-name=svim
#SBATCH --array=0-15
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/svim/logs/%A_%a.out
#SBATCH --error=/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/svim/logs/%A_%a.err

set -euo pipefail

# ── Set mode: "ont" or "pacbio" ──────────────────────────────────────────────
MODE="ont"   # <-- Change to "pacbio" for PacBio runs
             #     Also update --job-name above to "svim_pacbio"!
# ─────────────────────────────────────────────────────────────────────────────

source /usr/local/anaconda/5.1.0-Python3.6-gcc5/etc/profile.d/conda.sh
conda activate svim_env_clean

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
    OUTBASE="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/svim/svim_ont"
    bc_num="${barcode#bc}"
    bam="${BAM_DIR}/barcode${bc_num}_q10_mq20_dedup.bam"

elif [[ "${MODE}" == "pacbio" ]]; then
    BASE="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/demultiplex_pb/workflow_HiFiTargetEnrichment/batches/target_hifit"
    OUTBASE="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/svim/svim_pacbio"
    bam="${BASE}/${barcode}/realigned_minimap2/${barcode}.minimap2.GRCh38.bam"

else
    echo "ERROR: MODE must be 'ont' or 'pacbio'. Got: ${MODE}" >&2
    exit 1
fi

outdir="${OUTBASE}/${barcode}"

mkdir -p "${OUTBASE}/logs"
mkdir -p "${outdir}"

if [[ ! -f "${bam}" ]]; then
    echo "ERROR: BAM not found: ${bam}" >&2
    exit 1
fi

if [[ ! -f "${bam}.bai" ]]; then
    echo "ERROR: BAM index not found: ${bam}.bai" >&2
    exit 1
fi

echo "[$(date)] Running SVIM (${MODE}) for ${barcode} (task ${SLURM_ARRAY_TASK_ID})"
svim alignment \
    "${outdir}" \
    "${bam}" \
    "${REFERENCE}"
echo "[$(date)] Done: ${barcode}"