#!/bin/bash
#SBATCH --job-name=dedup_coordsorted_bam
#SBATCH --output=dedup_coordsorted_bam_%j.out
#SBATCH --error=dedup_coordsorted_bam_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

module load samtools

# --- Paths --- #
INPUT_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/bams_coord_sorted"
OUTPUT_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/bams_dedup_cleaned"

# --- Filters --- #
MIN_Q=10
MIN_MAPQ=20

# --- Create output directory --- #
mkdir -p "$OUTPUT_DIR"

# --- Summary file --- #
SUMMARY="${OUTPUT_DIR}/filtering_summary.tsv"
echo -e "barcode\tinput_reads\tQ_pass\tQ_MAPQ_pass\tfinal_dedup_reads" > "$SUMMARY"

# --- Define barcode list --- #
BARCODES=($(seq 1001 1008) $(seq 1097 1104))

# --- Loop through barcodes --- #
for barcode in "${BARCODES[@]}"; do
    echo "🔍 Processing barcode${barcode}..."

    IN_BAM="${INPUT_DIR}/barcode${barcode}_coordSorted.bam"
    FILTERED_BAM="${OUTPUT_DIR}/barcode${barcode}_q${MIN_Q}_mq${MIN_MAPQ}.bam"
    DEDUP_BAM="${OUTPUT_DIR}/barcode${barcode}_q${MIN_Q}_mq${MIN_MAPQ}_dedup.bam"

    if [[ ! -f "$IN_BAM" ]]; then
        echo "⚠️  File not found for barcode${barcode}, skipping..."
        continue
    fi

    # --- Read counts --- #
    INPUT_READS=$(samtools view -c "$IN_BAM")

    # Count Q-pass only
    Q_PASS=$(samtools view "$IN_BAM" \
        | awk -v minq="$MIN_Q" '
            {
                for (i=12; i<=NF; i++) {
                    if ($i ~ /^qs:f:/) {
                        split($i,a,":");
                        if (a[3] >= minq) print;
                        break
                    }
                }
            }' \
        | wc -l)

    # --- Apply Q + MAPQ filters in one pass --- #
    samtools view -h "$IN_BAM" \
    | awk -v minq="$MIN_Q" -v minmq="$MIN_MAPQ" '
        BEGIN { OFS="\t" }
        /^@/ { print; next }
        $5 < minmq { next }
        {
            for (i=12; i<=NF; i++) {
                if ($i ~ /^qs:f:/) {
                    split($i,a,":")
                    if (a[3] >= minq) print
                    break
                }
            }
        }' \
    | samtools view -b -o "$FILTERED_BAM"

    Q_MAPQ_PASS=$(samtools view -c "$FILTERED_BAM")

    # --- Deduplicate --- #
    samtools markdup -r "$FILTERED_BAM" "$DEDUP_BAM"
    samtools index "$DEDUP_BAM"

    FINAL_READS=$(samtools view -c "$DEDUP_BAM")

    # --- Save summary --- #
    echo -e "barcode${barcode}\t${INPUT_READS}\t${Q_PASS}\t${Q_MAPQ_PASS}\t${FINAL_READS}" >> "$SUMMARY"

    # --- Cleanup --- #
    rm -f "$FILTERED_BAM"

    echo "✅ Finished barcode${barcode}"
done

echo "🎯 Filtering (Q ≥ ${MIN_Q}, MAPQ ≥ ${MIN_MAPQ}) + deduplication complete!"
echo "📊 Summary written to ${SUMMARY}"
