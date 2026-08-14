#!/bin/bash
#SBATCH --job-name=bam2fastq
#SBATCH --output=bam2fastq_%j.out
#SBATCH --error=bam2fastq_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

module load samtools

BAM_DIR=/path/to/dorado_demultiplex_142pods_barcode_both_ends/bams_dedup_cleaned
FASTQ_DIR=/path/to/dorado_demultiplex_142pods_barcode_both_ends/fastq_dedup_cleaned

mkdir -p "${FASTQ_DIR}"

for bam in "${BAM_DIR}"/*.bam; do
    base=$(basename "${bam}" .bam)
    echo "Converting ${base}.bam → ${base}.fastq.gz"
    samtools fastq "${bam}" | gzip > "${FASTQ_DIR}/${base}.fastq.gz"
done

echo "All BAMs converted successfully."
