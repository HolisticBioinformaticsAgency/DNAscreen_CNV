#!/bin/bash
module load R

run_dnascreen_single() {
  local dnascreen_run=$1
  vcfs_and_bams_dir="/fs03/vh83/projects/temp_dnascreen_copy/dnascreen/consented/run${dnascreen_run}"
  sim_dir="/fs03/vh83/projects/temp_dnascreen_copy/dnascreen/exon_sims"
  sim_type="clinvar_sims"
  bed_name="9genes_25bp.fix.sorted"
  bed_file="${sim_dir}/bed_files/${bed_name}.bed"

  cd "$sim_dir/$sim_type" || exit
  if [ -d "run$dnascreen_run" ]; then
    echo "Directory run$dnascreen_run already exists."
  else
    mkdir "run$dnascreen_run"
    echo "Directory run$dnascreen_run created."
  fi

  if [ -d "run${dnascreen_run}/selected_regions" ]; then
    echo "Directory run${dnascreen_run}/selected_regions already exists."
  else
    mkdir "run$dnascreen_run/selected_regions"
    echo "Directory run$dnascreen_run/selected_regions created."
  fi

  # to avoid repeating 'test_dnascreen' R script when rerunning
  bams_rdata="${sim_dir}/RData/run${dnascreen_run}_bams_${bed_name}.RData"
  samples_with_cnvcalls="${sim_dir}/samples_with_cnvcalls/run${dnascreen_run}_samples_with_cnvcalls.txt"

  # Building vcf_lst.txt
  ls "${vcfs_and_bams_dir}/vcfs_run$dnascreen_run"/*.sorted.vcf.gz | xargs -n 1 basename > "run$dnascreen_run/vcf_lst.txt"

  # Check if bams_rdata exists
  if [ -f "$bams_rdata" ]; then
    echo "File $bams_rdata exists."
  else
    echo "File $bams_rdata does NOT exist."
  fi

  # Check if samples_with_cnvcalls exists
  if [ -f "$samples_with_cnvcalls" ]; then
    echo "File $samples_with_cnvcalls exists."
  else
    echo "File $samples_with_cnvcalls does NOT exist."
  fi

  # Make samples_with_cnvcalls.txt
  if [ ! -f "$bams_rdata" ] || [ ! -f "$samples_with_cnvcalls" ]; then
    Rscript "${sim_dir}/${sim_type}/master_script/test_dnascreen_run_different_ref_samples.R" \
    -dnascreen_run "$dnascreen_run" \
    -vcfs_and_bams_dir "$vcfs_and_bams_dir" \
    -sim_dir "$sim_dir" \
    -bed_file "$bed_file"
  else
    echo "The file $bams_rdata already exists. Skipping test_dnascreen_run_different_ref_samples.R."
  fi

  # Create samples_without_cnvcalls.txt
  vcf_lst="${sim_dir}/${sim_type}/run${dnascreen_run}/vcf_lst.txt"
  samples_without_cnvcalls="${sim_dir}/${sim_type}/run${dnascreen_run}/samples_without_cnvcalls.txt"

  # Extract base sample names (without .sorted.vcf.gz) for comparison
  awk -F'.' '{print $1}' "$vcf_lst" > "${sim_dir}/${sim_type}/run${dnascreen_run}/vcf_samples.txt"

  # Find the samples that are not in samples_with_cnvcalls.txt
  grep -F -x -v -f "$samples_with_cnvcalls" "${sim_dir}/${sim_type}/run${dnascreen_run}/vcf_samples.txt" | \
  while read -r sample; do
    # Directly use the sample name without the extension
    echo "$sample"
  done > "$samples_without_cnvcalls"

  # Check if samples_without_cnvcalls.txt was created successfully
  if [ -f "$samples_without_cnvcalls" ]; then
    echo "File $samples_without_cnvcalls created successfully."
  else
    echo "Failed to create $samples_without_cnvcalls."
  fi

  rm -rf "${sim_dir}/${sim_type}/run${dnascreen_run}/vcf_samples.txt"
}

export -f run_dnascreen_single


# Run for each dnascreen_run value in parallel
for dnascreen_run in 6 18 13 14 25 29; do
  run_dnascreen_single "$dnascreen_run" &
done

# Wait for all background processes to finish
wait

echo "All dnascreen_run jobs have completed."


# Function to run merge_sample_region_clinvar.sh for a pair of runs
run_merge_for_pair() {
  local run1=$1
  local run2=$2
  local sim_dir=$3
  local sim_type=$4

  echo "Processing pair: Run $run1 and Run $run2"

  # Call the merge_sample_region_clinvar.sh script for the pair of runs
  ${sim_dir}/${sim_type}/master_script/merge_sample_region_clinvar.sh "$sim_dir" "$run1" "$run2"
}

# Define the run pairs
run_pairs=(
  "6 18"
  "13 14"
  "25 29"
)


# Loop over the pairs
for pair in "${run_pairs[@]}"; do
  IFS=' ' read -ra run_pair <<< "$pair"
  run1="${run_pair[0]}"
  run2="${run_pair[1]}"
  sim_dir="/fs03/vh83/projects/temp_dnascreen_copy/dnascreen/exon_sims"
  sim_type="clinvar_sims"

  # Call the function to process the pair
  run_merge_for_pair "$run1" "$run2" "$sim_dir" "$sim_type"
done

echo "Merging samples and regions for all run pairs is completed."


run_simulations_for_run() {
  local dnascreen_run=$1
  sim_dir="/fs03/vh83/projects/temp_dnascreen_copy/dnascreen/exon_sims"
  vcfs_and_bams_dir="/fs03/vh83/projects/temp_dnascreen_copy/dnascreen/consented/run${dnascreen_run}"

  for cn in 1 3; do
    ${sim_dir}/${sim_type}/master_script/clinvar_simulation_repeat_cn${cn}.sh -dnascreen_run "$dnascreen_run" \
    -sim_dir "$sim_dir" -vcfs_and_bams_dir "$vcfs_and_bams_dir" &
  done

  Wait for both background jobs to finish
  wait

  if [ -d "${sim_dir}/clinvar_sims/run$dnascreen_run/plot_generation" ]; then
    echo "Directory run$dnascreen_run/plot_generation already exists."
  else
    mkdir -p "${sim_dir}/clinvar_sims/run$dnascreen_run/plot_generation/varseq_tables_highest_sensitivity/cov_samples"
    echo "Directory run$dnascreen_run/plot_generation created."
  fi

  # Make output_combined file for easy importing into varseq
  cd "${sim_dir}/clinvar_sims/run${dnascreen_run}" || exit
  mkdir -p output_combined
  cp output_cn1/* output_cn3/* output_combined/

  echo "All simulations for CN1 and CN3 for run $dnascreen_run have completed."
}

export -f run_simulations_for_run

# Run for each dnascreen_run value in parallel
for dnascreen_run in 6 18 13 14 25 29; do
  run_simulations_for_run "$dnascreen_run" &
done

# Wait for all background processes to finish
wait

echo "Simulations for all dnascreen_runs have completed."
