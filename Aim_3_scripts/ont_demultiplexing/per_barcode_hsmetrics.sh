#!/bin/bash
#SBATCH --job-name=hsmetrics_per_barcode
#SBATCH --output=hsmetrics_%j.out
#SBATCH --error=hsmetrics_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

module load picard

REF="/path/to/reference/genomes/hg38/heng_li_recomended/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
TARGETS="/path/to/dorado_demultiplex_142pods_barcode_both_ends/3539131_Covered_DNA_Screen_for_ONT_alignment.interval_list"
BAM_DIR="/path/to/dorado_demultiplex_142pods_barcode_both_ends/bams_demuxed"
OUT_DIR="/path/to/dorado_demultiplex_142pods_barcode_both_ends/hsmetrics"

mkdir -p "$OUT_DIR"

for BAM in "$BAM_DIR"/*barcode1*_rg.bam; do
    # Skip if file does not exist
    [[ ! -f "$BAM" ]] && continue

    SAMPLE=$(basename "$BAM" .bam)
    OUTFILE="${OUT_DIR}/${SAMPLE}.HSmetrics.txt"

    echo "Processing $SAMPLE ..."

    picard CollectHsMetrics \
        I="$BAM" \
        O="$OUTFILE" \
        R="$REF" \
        BAIT_INTERVALS="$TARGETS" \
        TARGET_INTERVALS="$TARGETS"
done

echo "✅ HsMetrics completed for all barcode BAMs."
