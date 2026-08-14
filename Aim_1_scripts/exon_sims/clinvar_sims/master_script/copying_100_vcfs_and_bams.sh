#!/bin/bash
# this is not runned anymore!

# Function to display usage information
usage() {
  echo "Usage: $0 -dnascreen_run <value> -exome_depth_dir \
  <value> -varseq_dir <value> -sim_dir <value>"
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
    -exome_depth_dir)
      exome_depth_dir="$2"
      shift # past argument
      shift # past value
      ;;
    -varseq_dir)
      varseq_dir="$2"
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
if [ -z "$dnascreen_run" ] || [ -z "$exome_depth_dir" ] || \
    [ -z "$varseq_dir" ] || [ -z "$sim_dir" ]; then
  usage
fi

original_vcf_dir="${varseq_dir}/vcfs_run${dnascreen_run}"
original_bam_dir="${exome_depth_dir}/bams_run${dnascreen_run}"

# Ensure the destination directory exists
mkdir -p "${sim_dir}/run${dnascreen_run}/vcfs_and_bams"

# Loop through each line in 100_vcf.txt
while IFS= read -r line; do
    # Extract the sample ID from the line
    sample_id=$(basename "$line" .sorted.vcf.gz)

    cp "${original_vcf_dir}/${line}" "${sim_dir}/run${dnascreen_run}/vcfs_and_bams/${line}"
    cp "${original_vcf_dir}/${line}.tbi" "${sim_dir}/run${dnascreen_run}/vcfs_and_bams/${line}.tbi"
    cp "${original_bam_dir}/${sample_id}.hq.sorted.marked.bam" "${sim_dir}/run${dnascreen_run}/vcfs_and_bams/${sample_id}.hq.sorted.marked.bam"
    cp "${original_bam_dir}/${sample_id}.hq.sorted.marked.bam.bai" "${sim_dir}/run${dnascreen_run}/vcfs_and_bams/${sample_id}.hq.sorted.marked.bam.bai"

done < "${exome_depth_dir}/simulation_sample_sets/run${dnascreen_run}_100_vcf.txt"