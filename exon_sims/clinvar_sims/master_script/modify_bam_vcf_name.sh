#!/bin/bash
module load samtools
module load bcftools

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <input_directory>"
  exit 1
fi

# Set input and output directories from command-line argument
INPUT_DIR="$1"
OUTPUT_DIR="${INPUT_DIR}_fixed_name"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Process each BAM file in the directory
for BAMFILE in "$INPUT_DIR"/*.bam; do
  # Extract the filename without the directory or extension
  FILENAME=$(basename "$BAMFILE" .bam)
  
  # Set the new sample name from the filename
  NEWSAMPLE="$FILENAME"
  
  # Create a temporary header file for the BAM
  samtools view -H "$BAMFILE" > header.txt
  
  # Modify the @RG line in the header file to match the new sample name
  sed -i "s/SM:[^[:space:]]*/SM:$NEWSAMPLE/" header.txt
  
  # Write the modified BAM file to the output directory
  samtools reheader header.txt "$BAMFILE" > "$OUTPUT_DIR/${FILENAME}.bam"
  
  # Remove the temporary header file
  rm header.txt
  
  # Generate the .bai index file for the new BAM in the output directory
  samtools index "$OUTPUT_DIR/${FILENAME}.bam"
  
  echo "Created updated BAM file: $OUTPUT_DIR/${FILENAME}.bam"
done

# Process each VCF file in the directory
for VCFFILE in "$INPUT_DIR"/*.vcf.gz; do
  # Extract the filename without the directory or extension
  FILENAME=$(basename "$VCFFILE" .vcf.gz)
  
  # Set the new sample name from the filename
  NEWSAMPLE="$FILENAME"
  
  # Use bcftools to replace the sample name in the VCF header and save to a new file in the output directory
  bcftools reheader -s <(echo "$NEWSAMPLE") -o "$OUTPUT_DIR/${FILENAME}.vcf.gz" "$VCFFILE"
  
  # Index the new VCF file in the output directory
  tabix -p vcf "$OUTPUT_DIR/${FILENAME}.vcf.gz"
  
  echo "Created updated VCF and index file: $OUTPUT_DIR/${FILENAME}.vcf.gz and $OUTPUT_DIR/${FILENAME}.vcf.gz.tbi"
done

echo "All BAM and VCF files processed and saved in $OUTPUT_DIR."
