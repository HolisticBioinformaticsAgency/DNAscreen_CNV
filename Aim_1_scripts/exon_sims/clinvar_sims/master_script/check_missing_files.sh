#!/bin/bash

# Define the folder path and required file extensions
folder_path="/Volumes/dnascreen_joshua/dnascreen/exon_sims/clinvar_sims/run18/output_combined_fixed_name"
required_extensions=("bam" "bam.bai" "sorted.vcf.gz" "sorted.vcf.gz.tbi")

# Extract unique sample IDs by stripping known extensions, ignoring '.bam.covtsf' files
sample_ids=($(ls "$folder_path" | grep -v '\.bam\.covtsf' | sed -E 's/\.(bam|bam\.bai|sorted\.vcf\.gz|sorted\.vcf\.gz\.tbi)$//' | sort -u))

# Check each sample ID for missing files
echo "Checking for missing files for each sample..."
for sample_id in "${sample_ids[@]}"; do
    missing_files=()

    # Check if each required extension exists for the current sample ID
    for ext in "${required_extensions[@]}"; do
        if [[ ! -f "$folder_path/${sample_id}.${ext}" ]]; then
            missing_files+=("$ext")
        fi
    done

    # Print missing files for each sample
    if (( ${#missing_files[@]} > 0 )); then
        echo "$sample_id is missing: ${missing_files[*]}"
    else
        echo "$sample_id has all required files."
    fi
done

# Print the final count of unique sample IDs
echo "Total unique sample IDs: ${#sample_ids[@]}"
