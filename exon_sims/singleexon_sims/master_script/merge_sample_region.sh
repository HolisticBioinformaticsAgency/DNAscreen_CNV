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

# Define the input files and directories
selected_region_dir="${sim_dir}/singleexon_sims/run${dnascreen_run}/selected_regions"
regions_files=("selected_regions_cn0" "selected_regions_cn1" "selected_regions_cn2" "selected_regions_cn3" "selected_regions_cn4")
samples_not_outlier_no_calls="${sim_dir}/singleexon_sims/run${dnascreen_run}/samples_not_outlier_no_calls.txt"

# Create an array to store sample IDs
mapfile -t sample_ids < "$samples_not_outlier_no_calls"

# Randomly select 166 samples (or the total number of samples if less than 166)
sample_count=${#sample_ids[@]}
if [ "$sample_count" -gt 166 ]; then
    selected_samples=($(shuf -n 166 -e "${sample_ids[@]}"))
else
    selected_samples=("${sample_ids[@]}")
fi

# Loop through each regions file
for regions_file in "${regions_files[@]}"
do
    # Create full path to regions file
    regions_file_path="${selected_region_dir}/${regions_file}"

    # Create output file path
    output_file="${selected_region_dir}/${regions_file}_sample"

    # Create a temporary file for storing modified output
    tmp_file="${output_file}.tmp"

    # Loop through both selected samples and bed_file line by line
    # Loop through both selected samples and bed_file line by line
    paste <(printf '%s\n' "${selected_samples[@]}" | sed 's/^[ \t]*//;s/[ \t]*$//') "$bed_file" | while IFS=$'\t' read -r sample_id bed_chr bed_start bed_end bed_gene; do
        # Trim leading/trailing spaces from the bed file columns
        bed_chr=$(echo "$bed_chr" | sed 's/^[ \t]*//;s/[ \t]*$//')
        bed_start=$(echo "$bed_start" | sed 's/^[ \t]*//;s/[ \t]*$//')
        bed_end=$(echo "$bed_end" | sed 's/^[ \t]*//;s/[ \t]*$//')
        bed_gene=$(echo "$bed_gene" | sed 's/^[ \t]*//;s/[ \t]*$//')

        # Create the new line with sample_id prefixed and the bed file coordinates
        new_line="${sample_id}\t${bed_chr}\t${bed_start}\t${bed_end}\t${bed_gene}"

        # Append the new line to the temporary file
        echo -e "$new_line" >> "$tmp_file"
    done


    # Overwrite original regions file with modified content
    mv "$tmp_file" "$output_file"

    # Display a message indicating completion for each regions file
    echo "Processed $regions_file_path. Output saved to $output_file"
done

echo "All files processed successfully."
