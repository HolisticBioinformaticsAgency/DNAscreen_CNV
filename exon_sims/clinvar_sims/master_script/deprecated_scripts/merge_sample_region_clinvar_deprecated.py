#!/usr/bin/env python3

import os
import sys

def read_samples(file_path):
    """Read samples from a file and return them as a list."""
    with open(file_path, 'r') as f:
        return [line.strip() for line in f]

def read_regions(file_path):
    """Read regions from a file and return them as a list."""
    with open(file_path, 'r') as f:
        return [line.strip() for line in f]

def write_output(output_path, samples, regions):
    """Write the samples and regions to the output file."""
    with open(output_path, 'w') as f:
        for i, region in enumerate(regions):
            sample = samples[i % len(samples)]
            f.write(f"{sample}\t{region}\n")

def process_run_pair(sim_dir, run1, run2):
    """Process a pair of runs and merge the sample-region data."""
    # Define file paths
    selected_regions_cn1 = os.path.join(sim_dir, 'clinvar_sims', 'clinvar_selected_regions', 'selected_regions_cn1')
    selected_regions_cn3 = os.path.join(sim_dir, 'clinvar_sims', 'clinvar_selected_regions', 'selected_regions_cn3')
    
    samples_run1 = os.path.join(sim_dir, 'clinvar_sims', f'run{run1}', 'samples_without_cnvcalls.txt')
    samples_run2 = os.path.join(sim_dir, 'clinvar_sims', f'run{run2}', 'samples_without_cnvcalls.txt')
    
    output_cn1_run1 = os.path.join(sim_dir, 'clinvar_sims', f'run{run1}', 'selected_regions', 'selected_regions_cn1_sample')
    output_cn1_run2 = os.path.join(sim_dir, 'clinvar_sims', f'run{run2}', 'selected_regions', 'selected_regions_cn1_sample')
    output_cn3_run1 = os.path.join(sim_dir, 'clinvar_sims', f'run{run1}', 'selected_regions', 'selected_regions_cn3_sample')
    output_cn3_run2 = os.path.join(sim_dir, 'clinvar_sims', f'run{run2}', 'selected_regions', 'selected_regions_cn3_sample')
    
    # Create output directories if they don't exist
    os.makedirs(os.path.dirname(output_cn1_run1), exist_ok=True)
    os.makedirs(os.path.dirname(output_cn1_run2), exist_ok=True)
    
    # Read samples and regions
    samples1 = read_samples(samples_run1)
    samples2 = read_samples(samples_run2)
    
    # Exit if no samples
    if not samples1 or not samples2:
        print(f"No samples found for runs {run1} or {run2}. Skipping pair.")
        return
    
    regions_cn1 = read_regions(selected_regions_cn1)
    regions_cn3 = read_regions(selected_regions_cn3)
    
    # Split regions for run1 and run2
    mid_index_cn1 = len(regions_cn1) // 2
    regions_cn1_run1 = regions_cn1[:mid_index_cn1]
    regions_cn1_run2 = regions_cn1[mid_index_cn1:]
    
    mid_index_cn3 = len(regions_cn3) // 2
    regions_cn3_run1 = regions_cn3[:mid_index_cn3]
    regions_cn3_run2 = regions_cn3[mid_index_cn3:]
    
    # Write output files
    write_output(output_cn1_run1, samples1, regions_cn1_run1)
    write_output(output_cn1_run2, samples2, regions_cn1_run2)
    write_output(output_cn3_run1, samples1, regions_cn3_run1)
    write_output(output_cn3_run2, samples2, regions_cn3_run2)
    
    print(f"Merged regions for run pair {run1} and {run2} into:")
    print(f"  {output_cn1_run1}")
    print(f"  {output_cn1_run2}")
    print(f"  {output_cn3_run1}")
    print(f"  {output_cn3_run2}")

def main():
    if len(sys.argv) != 4:
        print("Usage: merge_sample_region_clinvar.py <sim_dir> <run1> <run2>")
        sys.exit(1)

    # Parse command-line arguments
    sim_dir = sys.argv[1]
    run1 = sys.argv[2]
    run2 = sys.argv[3]
    
    # Process the run pair
    process_run_pair(sim_dir, run1, run2)

if __name__ == '__main__':
    main()
