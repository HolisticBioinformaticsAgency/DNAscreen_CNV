#!/bin/bash
#SBATCH --job-name=minimap2_pacbio
#SBATCH --output=minimap2_pacbio%j.out
#SBATCH --error=minimap2_pacbio%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G


set -euo pipefail


############################
# Load modules
############################
module load samtools/1.19.3
module load minimap2


############################
# SLURM resources
############################
THREADS=${SLURM_CPUS_PER_TASK}


############################
# Paths
############################
BASE_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/demultiplex_pb/workflow_HiFiTargetEnrichment/batches/target_hifit"
REF="/fs04/vh83/reference/genomes/hg38/heng_li_recomended/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"


############################
# QC thresholds
############################
MIN_MAPQ=20
MIN_READQUAL=10


############################
# Barcode list
############################
BARCODES=(
  bc1001 bc1002 bc1003 bc1004 bc1005 bc1006 bc1007 bc1008
  bc1097 bc1098 bc1099 bc1100 bc1101 bc1102 bc1103 bc1104
)


############################
# Loop over barcodes
############################
for BC in "${BARCODES[@]}"; do
  echo "===================================="
  echo "Processing ${BC}"
  echo "===================================="


  BC_DIR="${BASE_DIR}/${BC}"
  ALIGNED_BAM="${BC_DIR}/aligned/${BC}.GRCh38_noalt.bam"
  FASTQ_DIR="${BC_DIR}/fastq"
  REALIGN_DIR="${BC_DIR}/realigned_minimap2"


  mkdir -p "${FASTQ_DIR}" "${REALIGN_DIR}"


  FASTQ="${FASTQ_DIR}/${BC}.fastq"
  RAW_BAM="${REALIGN_DIR}/${BC}.minimap2.GRCh38.raw.bam"
  REALIGNED_BAM="${REALIGN_DIR}/${BC}.minimap2.GRCh38.bam"


  ############################
  # BAM -> FASTQ
  ############################
  echo "[${BC}] BAM → FASTQ"
  samtools fastq \
    -@ "${THREADS}" \
    "${ALIGNED_BAM}" > "${FASTQ}"


  ############################
  # Realign with minimap2
  ############################
  echo "[${BC}] minimap2 realignment"
  minimap2 \
    -t "${THREADS}" \
    -ax map-hifi \
    "${REF}" \
    "${FASTQ}" | \
    samtools sort \
      -@ "${THREADS}" \
      -o "${RAW_BAM}"


  ############################
  # Filter by mapping quality (MAPQ) and average read quality
  # Final output keeps the original naming convention
  ############################
  echo "[${BC}] Filtering: MAPQ ≥ ${MIN_MAPQ}, read quality ≥ ${MIN_READQUAL}"
  samtools view \
    -@ "${THREADS}" \
    -b \
    -q "${MIN_MAPQ}" \
    -e "avgqual(qual) >= ${MIN_READQUAL}" \
    -o "${REALIGNED_BAM}" \
    "${RAW_BAM}"


  samtools index "${REALIGNED_BAM}"


  ############################
  # Clean up intermediate raw BAM
  ############################
  rm -f "${RAW_BAM}"


  echo "[${BC}] Completed"
done


echo "All barcodes finished successfully 🎉"