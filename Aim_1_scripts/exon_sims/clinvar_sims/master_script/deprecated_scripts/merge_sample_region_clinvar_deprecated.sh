#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <sim_dir> <run1> <run2>"
  exit 1
fi

# Assign command-line arguments to variables
sim_dir="$1"
run1="$2"
run2="$3"

# Define input files
selected_regions_cn1="${sim_dir}/clinvar_sims/clinvar_selected_regions/selected_regions_cn1"
selected_regions_cn3="${sim_dir}/clinvar_sims/clinvar_selected_regions/selected_regions_cn3"

# Read samples from both runs
samples_run1="${sim_dir}/clinvar_sims/run${run1}/samples_without_cnvcalls.txt"
samples_run2="${sim_dir}/clinvar_sims/run${run2}/samples_without_cnvcalls.txt"

# Read samples into arrays
mapfile -t samples1 < "$samples_run1"
mapfile -t samples2 < "$samples_run2"

# Check if samples exist
if [ ${#samples1[@]} -eq 0 ] || [ ${#samples2[@]} -eq 0 ]; then
  echo "No samples found for runs $run1 or $run2. Skipping pair."
  exit 1
fi

# Initialize output files for run1 and run2
output_cn1_run1="${sim_dir}/clinvar_sims/run${run1}/selected_regions/selected_regions_cn1_sample"
output_cn1_run2="${sim_dir}/clinvar_sims/run${run2}/selected_regions/selected_regions_cn1_sample"
output_cn3_run1="${sim_dir}/clinvar_sims/run${run1}/selected_regions/selected_regions_cn3_sample"
output_cn3_run2="${sim_dir}/clinvar_sims/run${run2}/selected_regions/selected_regions_cn3_sample"

# Ensure the output directories exist
mkdir -p "${sim_dir}/clinvar_sims/run${run1}/selected_regions"
mkdir -p "${sim_dir}/clinvar_sims/run${run2}/selected_regions"

# Process CNVs for cn1
{
  # Read regions from selected_regions_cn1
  IFS=$'\n' read -d '' -r -a regions_cn1 < "$selected_regions_cn1"
  
  # Split the regions into two halves
  mid_index_cn1=$(( ${#regions_cn1[@]} / 2 ))
  regions_cn1_run1=("${regions_cn1[@]:0:mid_index_cn1}")  # First half for run1
  regions_cn1_run2=("${regions_cn1[@]:mid_index_cn1}")  # Second half for run2

  # Match samples with CNVs for run1
  for i in "${!regions_cn1_run1[@]}"; do
    sample="${samples1[i % ${#samples1[@]}]}"
    echo -e "$sample\t${regions_cn1_run1[$i]}"
  done
} > "$output_cn1_run1"

# Match samples with CNVs for run2
{
  for i in "${!regions_cn1_run2[@]}"; do
    sample="${samples2[i % ${#samples2[@]}]}"
    echo -e "$sample\t${regions_cn1_run2[$i]}"
  done
} > "$output_cn1_run2"

# Process CNVs for cn3
{
  # Read regions from selected_regions_cn3
  IFS=$'\n' read -d '' -r -a regions_cn3 < "$selected_regions_cn3"

  # Split the regions into two halves
  mid_index_cn3=$(( ${#regions_cn3[@]} / 2 ))
  regions_cn3_run1=("${regions_cn3[@]:0:mid_index_cn3}")  # First half for run1
  regions_cn3_run2=("${regions_cn3[@]:mid_index_cn3}")  # Second half for run2

  # Match samples with CNVs for run1
  for i in "${!regions_cn3_run1[@]}"; do
    sample="${samples1[i % ${#samples1[@]}]}"
    echo -e "$sample\t${regions_cn3_run1[$i]}"
  done
} > "$output_cn3_run1"

# Match samples with CNVs for run2
{
  for i in "${!regions_cn3_run2[@]}"; do
    sample="${samples2[i % ${#samples2[@]}]}"
    echo -e "$sample\t${regions_cn3_run2[$i]}"
  done
} > "$output_cn3_run2"

echo "Merged regions for run pair $run1 and $run2 into:"
echo "$output_cn1_run1"
echo "$output_cn1_run2"
echo "$output_cn3_run1"
echo "$output_cn3_run2"
