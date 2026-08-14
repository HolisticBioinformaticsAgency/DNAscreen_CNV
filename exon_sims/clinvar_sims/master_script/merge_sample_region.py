import argparse
import pandas as pd
import os
import numpy as np

# Function to display usage information
def usage():
    print("Usage: script.py -dnascreen_run <value> -sim_dir <value> -sim_type <value> -bed_file <value> -padding <value> -cn <values>")
    exit(1)

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Merge samples with BED file.")
    parser.add_argument('-dnascreen_run', required=True, help='DNAscreen run value')
    parser.add_argument('-sim_dir', required=True, help='Simulation directory')
    parser.add_argument('-sim_type', required=True, help='Simulation type')
    parser.add_argument('-bed_file', required=True, help='Path to the BED file')
    parser.add_argument('-padding', type=int, default=0, help='Padding size for BED regions (default: 75)')
    parser.add_argument('-cn', nargs='+', type=int, default=[0, 1, 2, 3, 4],
                        help='List of CN values to pick (default: 0 1 2 3 4)')

    args = parser.parse_args()

    # Define the input files and directories
    sim_dir = args.sim_dir
    sim_type = args.sim_type
    dnascreen_run = args.dnascreen_run
    padding = args.padding
    selected_region_dir = os.path.join(sim_dir, f"{sim_type}/run{dnascreen_run}/selected_regions")
    
    # List of regions files based on selected CN values
    regions_files = [f"selected_regions_cn{cn}" for cn in args.cn]
    
    samples_file_path = os.path.join(sim_dir, f"{sim_type}/run{dnascreen_run}/samples_not_outlier_no_calls.txt")
    
    # Read the samples_not_outlier_no_calls file into a DataFrame
    samples_df = pd.read_csv(samples_file_path, header=None, names=['sample_id'], squeeze=True)

    # Read the BED file into a DataFrame
    bed_df = pd.read_csv(args.bed_file, delim_whitespace=True, header=None)

    # Loop through each regions file specified by the user
    for regions_file in regions_files:
        # Create output file path
        output_file = os.path.join(selected_region_dir, f"{regions_file}_sample")

        # Initialize a list to hold the new data
        new_data = []

        # Ensure random sampling without replacement
        np.random.seed(40)  # Optional: set seed for reproducibility # previous is 42
        sample_ids = samples_df.tolist()

        # Loop through each row in the BED file and apply padding
        for _, row in bed_df.iterrows():
            # Check if the row has at least 4 columns
            if len(row) < 4:
                print(f"Warning: Row has fewer than 4 columns. Skipping row: {row}")
                continue

            # Extracting bed file columns and convert to string, adding padding
            bed_chr = str(row[0]).strip()
            bed_start = max(0, int(row[1]) - padding)  # Ensure start is not negative
            bed_end = int(row[2]) + padding
            bed_gene = str(row[3]).strip()

            # Randomly pick a sample from the samples_df without replacement
            if sample_ids:  # Check if there are still samples available
                sample_id = np.random.choice(sample_ids, replace=False)
                sample_ids.remove(sample_id)  # Remove the chosen sample to avoid replacement

                # Create the new line with sample_id prefixed and the bed file coordinates
                new_line = f"{sample_id}\t{bed_chr}\t{bed_start}\t{bed_end}\t{bed_gene}"
                new_data.append(new_line)
            else:
                print("Warning: No more samples available for selection.")
                break

        # Write the new data to the output file, with an additional newline at the end
        with open(output_file, 'w') as f:
            f.write("\n".join(new_data))
            f.write("\n")  # Add an empty line at the bottom of the file

        print(f"Processed {regions_file}. Output saved to {output_file}")

if __name__ == "__main__":
    main()
