#!/bin/bash
#SBATCH --job-name=svim
#SBATCH --array=0-15
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/path/to/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/svim/logs/%A_%a.out
#SBATCH --error=/path/to/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/svim/logs/%A_%a.err


set -euo pipefail


# ── Set mode: "ont" or "pacbio" ──────────────────────────────────────────────
MODE="ont"   # <-- Change to "pacbio" for PacBio runs
             #     Also update --job-name above to "svim_pacbio"!
# ─────────────────────────────────────────────────────────────────────────────


source /usr/local/anaconda/5.1.0-Python3.6-gcc5/etc/profile.d/conda.sh
conda activate svim_env_clean


REFERENCE="/path/to/reference/genomes/hg38/heng_li_recomended/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"


# Full list of round2 barcodes (indices 0-15): bc1113-bc1128
BARCODES=(
    bc1113 bc1114 bc1115 bc1116 bc1117 bc1118 bc1119 bc1120
    bc1121 bc1122 bc1123 bc1124 bc1125 bc1126 bc1127 bc1128
)


barcode=${BARCODES[$SLURM_ARRAY_TASK_ID]}
THREADS=${SLURM_CPUS_PER_TASK}


if [[ "${MODE}" == "ont" ]]; then
    BAM_DIR="/path/to/dorado_demultiplex_round2_full_ont/bams_dedup_cleaned"
    OUTBASE="/path/to/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/svim/svim_ont_round2"
    bc_num="${barcode#bc}"
    bam="${BAM_DIR}/barcode${bc_num}_q10_mq20_dedup.bam"


elif [[ "${MODE}" == "pacbio" ]]; then
    BASE="/path/to/demultiplex_pb/workflow_HiFiTargetEnrichment/batches/target_hifit"
    OUTBASE="/path/to/ONT_PacBio_CNV_calling/ont_pacbio_sv_callers/svim/svim_pacbio"
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
