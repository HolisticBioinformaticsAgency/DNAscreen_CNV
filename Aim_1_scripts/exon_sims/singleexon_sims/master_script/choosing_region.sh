#!/bin/bash

# Function to display usage information
usage() {
  echo "Usage: $0 -dnascreen_run <value> -sim_dir <value> -bed_file <value>"
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
      -bed_file)
      bed_file="$2"
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
if [ -z "$dnascreen_run" ] || [ -z "$sim_dir" ] || [ -z "$bed_file" ]; then
  usage
fi

# Define the output directory
output_dir="${sim_dir}/singleexon_sims/run${dnascreen_run}/selected_regions"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Define the output file names
output_files=("selected_regions_cn0" "selected_regions_cn1" "selected_regions_cn2" "selected_regions_cn3" "selected_regions_cn4")

# Loop through each output file
for output_file in "${output_files[@]}"; do
    # Overwrite the file by redirecting output to it
    > "$output_dir/$output_file"

    # Generate 20 CNV intervals
    for ((i = 1; i <= 20; i++)); do
    
      # Randomly select a line from the bed file
      random_bed_line=$(shuf -n 1 "$bed_file")
      echo "Randomly selected region: $random_bed_line"

      # Extract chromosome, start, and end positions from random_bed_line
      chrom=$(echo "$random_bed_line" | awk '{print $1}')
      start=$(echo "$random_bed_line" | awk '{print $2}')
      end=$(echo "$random_bed_line" | awk '{print $3}')
      reg_info=$(echo "$random_bed_line" | awk '{print $4}')

      # Adjust start and end positions by padding with +- 74 bases
      # i did this because reads that aren't fully within the window will not be
      # down-sampled or removed by bedtools
      adjusted_start=$((start - 74))
      adjusted_end=$((end + 74))

      # Construct the adjusted region string
      region="$chrom\t$adjusted_start\t$adjusted_end\t$reg_info"

      # # Construct the adjusted region string
      # region="$chrom\t$start\t$end\t$reg_info"

      # Output the chosen exon to the current output file
      echo -e "$region" >> "$output_dir/$output_file"
    done
done
