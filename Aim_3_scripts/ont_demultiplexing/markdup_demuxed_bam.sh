#!/bin/bash
#SBATCH --job-name=markdup_demuxed_bam
#SBATCH --output=markdup_demuxed_bam_%j.out
#SBATCH --error=markdup_demuxed_bam_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

module load samtools

# Paths
INPUT_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/bams_demuxed"
MARKDUP_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/bams_markdup"
COORDSORTED_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/bams_coord_sorted"
BED_FILE="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/3539131_Covered_DNA_Screen_for_ONT_alignment.bed"
BED_BAM_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/bams_bed_dan_restricted"

# Create output directories
mkdir -p "$MARKDUP_DIR"
mkdir -p "$COORDSORTED_DIR"
mkdir -p "$BED_BAM_DIR"

# Output TSV file
OUT_TSV="$MARKDUP_DIR/duplication_rates.tsv"
echo -e "Barcode\tDuplicationRate" > "$OUT_TSV"

# --- Define barcode ranges --- #
BARCODES=($(seq 1001 1008) $(seq 1097 1104))

# --- Loop through barcodes --- #
for barcode in "${BARCODES[@]}"; do
    echo "Processing barcode$barcode..."

    INPUT_BAM="$INPUT_DIR/PBG10946_pass_fb6074ec_f8186473_0_custom_barcode_0_19_custom_barcode_barcode${barcode}_rg.bam"
    BED_BAM="$BED_BAM_DIR/barcode${barcode}_bedRestricted.bam"
    COORDSORT_BAM="$COORDSORTED_DIR/barcode${barcode}_coordSorted.bam"
    MARKDUP_BAM="$MARKDUP_DIR/markdup_barcode${barcode}.bam"

    # 1. Sort by coordinate
    samtools sort "$INPUT_BAM" -o "$COORDSORT_BAM"
    samtools index "$COORDSORT_BAM"

    # 2. Restrict to BED regions
    samtools view -b -L "$BED_FILE" "$COORDSORT_BAM" > "$BED_BAM"

    # 3. Mark duplicates
    samtools markdup "$BED_BAM" "$MARKDUP_BAM"

    # 4. Calculate duplication rate
    dup_rate=$(echo "scale=4; $(samtools view -c -f 1024 $MARKDUP_BAM) / $(samtools view -c -F 4 $MARKDUP_BAM)" | bc)

    # 5. Append to TSV
    echo -e "barcode${barcode}\t$dup_rate" >> "$OUT_TSV"
done

echo "Done! Duplication rates saved to $OUT_TSV"
