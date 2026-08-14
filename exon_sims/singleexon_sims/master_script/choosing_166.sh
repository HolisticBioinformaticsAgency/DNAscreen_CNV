#!/bin/bash

# Function to display usage information
usage() {
  echo "Usage: $0 -dnascreen_run <value> -sim_dir <value> -vcfs_and_bams_dir <value>"
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
    *)
      echo "Unknown option: $key"
      usage
      ;;
  esac
done

# Check if all required arguments are provided
if [ -z "$dnascreen_run" ] || [ -z "$sim_dir" ]; then
  usage
fi

# Input file names
samples_not_outlier_no_calls="${sim_dir}/singleexon_sims/run${dnascreen_run}/samples_not_outlier_no_calls.txt"
vcf_list="${sim_dir}/singleexon_sims/run${dnascreen_run}/vcf_lst.txt"

# Output file name
output_file="${sim_dir}/singleexon_sims/run${dnascreen_run}/selected_166_vcf.txt"

# Check if there are at least 166 samples in the input file
sample_count=$(wc -l < "$samples_not_outlier_no_calls")
if [ "$sample_count" -lt 166 ]; then
  echo "Not enough samples. Found: $sample_count, Required: 166"
  exit 1
fi

# Create an array to hold the selected VCF files
selected_vcfs=()

# Loop until we have 166 valid VCF files
while [ ${#selected_vcfs[@]} -lt 166 ]; do
  # Randomly select a sample from the samples file
  sample=$(shuf -n 1 "$samples_not_outlier_no_calls")

  # Construct the expected VCF filename based on the sample ID
  vcf_file="${sample}.sorted.vcf.gz"

  # Check if the VCF file exists in the vcf_list
  if grep -q "$(basename "$vcf_file")" "$vcf_list"; then
    echo "$vcf_file exists!"
    selected_vcfs+=("$vcf_file")
  fi
done

# Write selected VCF files to the output file, one per line
printf "%s\n" "${selected_vcfs[@]}" > "$output_file"

# Confirm the output file has exactly 166 lines
final_line_count=$(wc -l < "$output_file")
if [ "$final_line_count" -eq 166 ]; then
  echo "Successfully selected 166 VCF files written to $output_file"
else
  echo "Error: Expected 166 VCF files, but found $final_line_count"
  exit 1
fi
