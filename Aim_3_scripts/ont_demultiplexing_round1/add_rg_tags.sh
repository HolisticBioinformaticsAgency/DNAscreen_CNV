#!/bin/bash
# Add missing @RG SM tags to BAM files based on specific barcode list
# Dependencies: samtools

module load samtools

# Directory containing BAM files
BAM_DIR="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/bams_demuxed"
cd "$BAM_DIR" || exit

# Barcode list
BARCODES=($(seq 1001 1008) $(seq 1097 1104))

for barcode in "${BARCODES[@]}"; do
    # Find the BAM file corresponding to this barcode
    bam_file=$(ls *barcode${barcode}*.bam 2>/dev/null)
    
    if [[ -z "$bam_file" ]]; then
        echo "⚠️  No BAM found for barcode${barcode}, skipping..."
        continue
    fi

    # Define output BAM filename with _rg suffix
    out_bam="${bam_file%.bam}_rg.bam"

    echo "Processing $bam_file ..."

    # Add @RG line (SM tag = barcode number)
    samtools addreplacerg \
        -r "@RG\tID:${barcode}\tSM:${barcode}\tPL:ONT\tLB:DNA_Screen" \
        -o "$out_bam" \
        "$bam_file"

    # Index new BAM
    samtools index "$out_bam"

    echo "✅ Added RG tag to $bam_file → $out_bam"
done

echo "All specified BAM files processed successfully."
