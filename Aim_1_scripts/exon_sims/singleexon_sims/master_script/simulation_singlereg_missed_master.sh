#!/bin/bash
module load R

# Specify a single run or a list of runs here
# Example: dnascreen_runs=(33 34 35)
dnascreen_runs=(13 14 29)

vcfs_and_bams_dir_base="/fs03/vh83/projects/temp_dnascreen_copy/dnascreen/consented"
sim_dir="/fs03/vh83/projects/temp_dnascreen_copy/dnascreen/exon_sims"
sim_type="singleexon_sims"
padding=0
prev_run_folder_name="_missed"
work_folder_name="_missed_2" #default is ""

run_dnascreen() {
  local dnascreen_run=$1
  local vcfs_and_bams_dir="${vcfs_and_bams_dir_base}/run${dnascreen_run}"
  missed_cnvs="${sim_dir}/${sim_type}/run${dnascreen_run}${prev_run_folder_name}/run${dnascreen_run}${work_folder_name}_singleexon.txt"

  # Read cn values from missed_cnvs file (extract unique CN values)
  cn_values=($(awk '{print $1}' "$missed_cnvs" | sort -u))

  # # Create the BED file (exclude the CN column)
  # bed_file="${sim_dir}/bed_files/run${dnascreen_run}${work_folder_name}_singleexon.bed"
  # awk '{print $2"\t"$3"\t"$4"\t"$5}' "$missed_cnvs" > "$bed_file"

  # Create necessary directories
  cd "$sim_dir/${sim_type}" || exit
  if [ -d "run${dnascreen_run}${work_folder_name}" ]; then
    echo "Directory run${dnascreen_run}${work_folder_name} already exists."
  else
    mkdir "run${dnascreen_run}${work_folder_name}"
    echo "Directory run${dnascreen_run}${work_folder_name} created."
  fi

  # Check if the samples_not_outlier_no_calls.txt file exists and copy it to the current work folder
  if [ -f "${sim_dir}/${sim_type}/run${dnascreen_run}/samples_not_outlier_no_calls.txt" ]; then
    echo "Copying samples_not_outlier_no_calls.txt to current work folder."
    cp "${sim_dir}/${sim_type}/run${dnascreen_run}/samples_not_outlier_no_calls.txt" \
       "${sim_dir}/${sim_type}/run${dnascreen_run}${work_folder_name}/samples_not_outlier_no_calls.txt"
  else
    echo "samples_not_outlier_no_calls.txt does not exist in run${dnascreen_run}."
  fi

  # Path to samples_not_outlier_no_calls
  samples_not_outlier_no_calls="${sim_dir}/${sim_type}/run${dnascreen_run}/samples_not_outlier_no_calls.txt"

  # Building vcf_lst.txt
  ls "${vcfs_and_bams_dir}/vcfs_run$dnascreen_run"/*.sorted.vcf.gz | xargs -n 1 basename > "run$dnascreen_run/vcf_lst.txt"

  # Check if samples_not_outlier_no_calls exists
  if [ -f "$samples_not_outlier_no_calls" ]; then
    echo "File $samples_not_outlier_no_calls exists."
  else
    echo "File $samples_not_outlier_no_calls does NOT exist."
  fi

  # Create folder for selected regions
  mkdir -p "${sim_dir}/${sim_type}/run${dnascreen_run}${work_folder_name}/selected_regions"

  # Run Python script for matching samples to regions
  for cn in "${cn_values[@]}"; do
    # Filter regions for the current CN value
    filtered_bed_file="${sim_dir}/${sim_type}/run${dnascreen_run}${work_folder_name}/selected_regions/bed_file_cn${cn}.bed"
    awk -v cn="$cn" '$1 == cn {print $2"\t"$3"\t"$4"\t"$5}' "$missed_cnvs" > "$filtered_bed_file"

    # Run the merge_sample_region.py script for the current CN
    python ./master_script/merge_sample_region.py \
    -dnascreen_run "${dnascreen_run}${work_folder_name}" \
    -sim_dir "$sim_dir" \
    -sim_type "$sim_type" \
    -bed_file "$filtered_bed_file" \
    -padding "$padding" \
    -cn "$cn"
  done

  # Run simulation_repeat_cn.sh scripts in parallel
  for cn in "${cn_values[@]}"; do
    ./master_script/simulation_repeat_cn${cn}.sh \
    -dnascreen_run "$dnascreen_run" \
    -work_folder_name "$work_folder_name" \
    -sim_dir "$sim_dir" \
    -sim_type "$sim_type" \
    -vcfs_and_bams_dir "$vcfs_and_bams_dir" &
  done

  # Wait for all background jobs to finish
  wait

  echo "All simulation_repeat_cn.sh scripts have completed."

  if [ -d "run${dnascreen_run}${work_folder_name}/plot_generation" ]; then
    echo "Directory run$dnascreen_run/plot_generation already exists."
  else
    mkdir -p "run${dnascreen_run}${work_folder_name}/plot_generation/varseq_tables/cov_samples"
    echo "Directory run$dnascreen_run/plot_generation created."
  fi

  # Combine outputs for VarSeq
  cd "run${dnascreen_run}${work_folder_name}" || exit
  mkdir -p output_combined

  # Loop over each CN value and copy files to output_combined
  for cn in "${cn_values[@]}"; do
    cp "output_cn${cn}/"* output_combined/
  done

  # Modify BAM and VCF headers for VarSeq compatibility
  ${sim_dir}/${sim_type}/master_script/modify_bam_vcf_name.sh "${sim_dir}/${sim_type}/run${dnascreen_run}${work_folder_name}/output_combined"

  # Return to simulation directory
  cd "$sim_dir/${sim_type}" || exit
}

export -f run_dnascreen

# Run simulations for each dnascreen_run in parallel
for dnascreen_run in "${dnascreen_runs[@]}"; do
  run_dnascreen "$dnascreen_run" &
done

# Wait for all background jobs to finish
wait

echo "All dnascreen_run jobs have completed."

# File count validation loop
echo "Starting file count validation..."
for dnascreen_run in "${dnascreen_runs[@]}"; do
  missed_cnvs="${sim_dir}/bed_files/run${dnascreen_run}${work_folder_name}_singleexon.txt"
  # Number of intervals in the bed file
  num_cnvs=$(wc -l "$missed_cnvs")
  expected_file_count=$num_cnvs
  actual_file_count=$(find "${sim_dir}/${sim_type}/run${dnascreen_run}${work_folder_name}/output_combined_fixed_name" -type f | wc -l)

  if [ "$actual_file_count" -eq "$expected_file_count" ]; then
    echo "Run ${dnascreen_run}: The output_combined directory contains the correct number of files: $actual_file_count files."
  else
    echo "Run ${dnascreen_run}: Warning - Incorrect number of files in output_combined."
    echo "Expected: $expected_file_count files, but found: $actual_file_count files."
  fi
done

echo "File count validation completed."
