#!/bin/bash

# Set input and output directories
INPUT_DIR="/Volumes/dnascreen_joshua/dnascreen/exon_sims/clinvar_sims/run6/output_combined"
OUTPUT_DIR="/Volumes/dnascreen_joshua/dnascreen/exon_sims/clinvar_sims/run6/output_combined_fixed_name"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each VCF file in the directory
for VCFFILE in "$INPUT_DIR"/*.vcf.gz; do
  # Extract the filename without the directory or extension
  FILENAME=$(basename "$VCFFILE" .vcf.gz)
  
  # Set the new sample name from the filename (including 'cn1_')
  NEWSAMPLE="$FILENAME"
  
  # Use bcftools to replace the sample name in the VCF header and save to a new file
  bcftools reheader -s <(echo "$NEWSAMPLE") -o "$OUTPUT_DIR/${FILENAME}.vcf.gz" "$VCFFILE"
  
  # Index the new VCF file and create a new .tbi file in the output directory
  tabix -p vcf "$OUTPUT_DIR/${FILENAME}.vcf.gz"
  
  echo "Created new VCF and index file: $OUTPUT_DIR/${FILENAME}.vcf.gz and $OUTPUT_DIR/${FILENAME}.vcf.gz.tbi"
done

echo "All VCF files processed and saved in $OUTPUT_DIR."
