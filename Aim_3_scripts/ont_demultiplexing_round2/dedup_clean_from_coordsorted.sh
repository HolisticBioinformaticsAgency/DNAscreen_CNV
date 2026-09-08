#!/bin/bash
#SBATCH --job-name=dedup_coordsorted_bam_round2_full
#SBATCH --output=dedup_coordsorted_bam_round2_full_%j.out
#SBATCH --error=dedup_coordsorted_bam_round2_full_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

module load samtools

# --- Paths --- #
# Input BAMs are already renamed + coordinate-sorted by rename_part.sh
INPUT_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_round2_full_ont/bams_demuxed"
OUTPUT_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_round2_full_ont/bams_dedup_cleaned"

# --- Filters --- #
MIN_Q=10
MIN_MAPQ=20

# --- Create output directory --- #
mkdir -p "$OUTPUT_DIR"

# --- Summary file --- #
SUMMARY="${OUTPUT_DIR}/filtering_summary.tsv"
echo -e "barcode\tinput_reads\tQ_pass\tQ_MAPQ_pass\tfinal_dedup_reads" > "$SUMMARY"

# --- Define barcode list --- #
# BCF/BCR barcodes 1113-1128
BARCODES=($(seq 1113 1128))

# --- Loop through barcodes --- #
for barcode in "${BARCODES[@]}"; do
    echo "🔍 Processing barcode${barcode}..."

    IN_BAM="${INPUT_DIR}/PBG10946_pass_fb6074ec_f8186473_0_custom_barcode_barcode${barcode}_rg.bam"
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
