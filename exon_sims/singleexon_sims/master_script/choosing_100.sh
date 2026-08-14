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
    -vcfs_and_bams_dir)
      vcfs_and_bams_dir="$2"
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
if [ -z "$dnascreen_run" ] || [ -z "$sim_dir" ] || [ -z "$vcfs_and_bams_dir" ]; then
  usage
fi

# Input file names
samples_with_cnvcalls="${sim_dir}/samples_with_cnvcalls/run${dnascreen_run}_samples_with_cnvcalls.txt"
vcf_list="${sim_dir}/singleexon_sims/run${dnascreen_run}/vcf_lst.txt"

# Output file name
output_file="${sim_dir}/simulation_sample_sets/run${dnascreen_run}_100_vcf.txt"

# Get a list of all sample IDs from samples_with_cnvcalls
sample_ids=$(awk '{print $1}' "$samples_with_cnvcalls")

# Get a list of all VCF files in vcf_lst.txt (reading line by line to handle spaces)
filtered_vcf_files=()
while IFS= read -r vcf_file; do
    # Extract the sample ID from the VCF filename (assuming filename format is 'sample_id.sorted.vcf.gz')
    sample_id=$(echo "$vcf_file" | awk -F'.' '{print $1}')
    
    # Check if the extracted sample ID is NOT in the samples_with_cnvcalls
    if ! echo "$sample_ids" | grep -q "$sample_id"; then
        filtered_vcf_files+=("$vcf_file")
    fi
done < "$vcf_list"

# Count the number of available VCF files after filtering
num_vcf_files=${#filtered_vcf_files[@]}

# If there are fewer than 100 valid VCF files, display a message and exit
if [ "$num_vcf_files" -lt 100 ]; then
    echo "Error: Less than 100 VCF files available after filtering."
    exit 1
fi

# Randomly select 100 VCF files from the filtered list
selected_vcf_files=()
for ((i=0; i<100; i++)); do
    # Choose a random index within the range of available files
    random_index=$((RANDOM % num_vcf_files))
    
    # Get the VCF filename at the random index
    selected_vcf_file="${filtered_vcf_files[random_index]}"
    
    # Append the selected filename to the list
    selected_vcf_files+=("$selected_vcf_file")
    
    # Remove the selected filename from the list to avoid duplicate selection
    unset 'filtered_vcf_files[random_index]'
    
    # Reset array to remove null entries
    filtered_vcf_files=("${filtered_vcf_files[@]}")
    
    # Update the count of available VCF files
    num_vcf_files=${#filtered_vcf_files[@]}
done

# Write the selected VCF files to the output file (each on a new line)
printf "%s\n" "${selected_vcf_files[@]}" > "$output_file"

echo "Randomly selected 100 VCF files that do not have sample IDs listed in $samples_with_cnvcalls."
echo "Output written to $output_file."
