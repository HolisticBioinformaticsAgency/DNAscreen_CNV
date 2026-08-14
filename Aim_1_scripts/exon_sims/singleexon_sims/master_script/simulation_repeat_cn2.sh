#!/bin/bash
module load bedtools
module load samtools
module load picard

# Function to display usage information
usage() {
  echo "Usage: $0 -dnascreen_run <value> -sim_dir <value> -sim_type <value> -vcfs_and_bams_dir <value> -work_folder_name (optional) <value>"
  exit 1
}

# Parse named arguments using getopts
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -dnascreen_run)
      dnascreen_run="$2"
      shift # past argument
      shift # past value
      ;;
    -sim_dir)
      sim_dir="$2"
      shift # past argument
      shift # past value
      ;;
    -sim_type)
      sim_type="$2"
      shift # past argument
      shift # past value
      ;;
    -vcfs_and_bams_dir)
      vcfs_and_bams_dir="$2"
      shift # past argument
      shift # past value
      ;;
    -work_folder_name)
      work_folder_name="$2"
      shift # past argument
      shift # past value
      ;;
    *)
      echo "Unknown option: $key"
      usage
      ;;
  esac
done

# Check if all required arguments are provided
if [ -z "$dnascreen_run" ] || [ -z "$sim_dir" ] || [ -z "$sim_type" ] || [ -z "$vcfs_and_bams_dir" ]; then
  usage
fi

# Allow empty work_folder_name and handle it
if [ -z "$work_folder_name" ]; then
  echo "No work folder name provided, proceeding without it."
  work_folder_name=""
fi

# Proceed with the rest of the script
echo "Working with folder extension: $work_folder_name"

# Set up directories and file paths
cnstate="2"
bams_dir="${vcfs_and_bams_dir}/bams_run${dnascreen_run}"
vcfs_dir="${vcfs_and_bams_dir}/vcfs_run${dnascreen_run}"
output_dir="${sim_dir}/${sim_type}/run${dnascreen_run}${work_folder_name}/output_cn${cnstate}"
selected_regions="${sim_dir}/${sim_type}/run${dnascreen_run}${work_folder_name}/selected_regions/selected_regions_cn${cnstate}_sample"

# Create output directory if it does not exist
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

# Loop through each line in selected_regions file
while IFS=$'\t' read -r sample_id chrom start end gene; do
    # Skip header lines or empty lines
    if [[ "$sample_id" == "Sample_ID" || -z "$sample_id" ]]; then
        continue
    fi

    echo "Processing sample: $sample_id"

    # Find corresponding VCF file
    sample_vcf_file="$vcfs_dir/${sample_id}.sorted.vcf.gz"
    if [ ! -f "$sample_vcf_file" ]; then
        echo "VCF file not found for sample: $sample_id"
        continue
    fi

    # Copy corresponding VCF file to output directory
    cp "$sample_vcf_file" "$output_dir/cn${cnstate}_${sample_id}.sorted.vcf.gz"
    cp "${sample_vcf_file}.tbi" "$output_dir/cn${cnstate}_${sample_id}.sorted.vcf.gz.tbi"

    # Find corresponding BAM file
    sample_bam_file="$bams_dir/${sample_id}.hq.sorted.marked.bam"
    if [ ! -f "$sample_bam_file" ]; then
        echo "BAM file not found for sample: $sample_id"
        continue
    fi

    # Window of CNV for upsampling
    window="${chrom}	${start}	${end}"
    echo "Simulated CNV window: $window"

    # Extract reads from the selected region (removed -wa option 6/11/24)
    bedtools intersect -abam "$sample_bam_file" -b <(echo "$window") > "$output_dir/extracted_$sample_id.bam"
    samtools index "$output_dir/extracted_$sample_id.bam"

    # Create original.bam excluding the selected region
    bedtools intersect -abam "$sample_bam_file" -b <(echo "$window") -v > "$output_dir/original_$sample_id.bam"
    samtools index "$output_dir/original_$sample_id.bam"

    # Downsample extracted.bam (P=1.0 retaining probability) 
    samtools view -b "$output_dir/extracted_$sample_id.bam" > "$output_dir/downsampled_$sample_id.bam"
    samtools index "$output_dir/downsampled_$sample_id.bam"

    # Merge downsampled.bam with original.bam
    picard MergeSamFiles INPUT="$output_dir/downsampled_$sample_id.bam" INPUT="$output_dir/original_$sample_id.bam" OUTPUT="$output_dir/cn${cnstate}_$sample_id.bam"
    samtools index "$output_dir/cn${cnstate}_$sample_id.bam"

    rm -rf "$output_dir/extracted_$sample_id.bam"
    rm -rf "$output_dir/extracted_$sample_id.bam.bai"
    rm -rf "$output_dir/original_$sample_id.bam"
    rm -rf "$output_dir/original_$sample_id.bam.bai"
    rm -rf "$output_dir/downsampled_$sample_id.bam"
    rm -rf "$output_dir/downsampled_$sample_id.bam.bai"


done < "$selected_regions"

echo "Script completed successfully."
