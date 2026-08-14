library(ExomeDepth)
library(tools)

# Function to parse named command line arguments
parse_args <- function(args) {
  parsed_args <- list()
  for (i in seq(1, length(args), by = 2)) {
    key <- gsub("^-", "", args[i])
    value <- args[i + 1]
    parsed_args[[key]] <- value
  }
  return(parsed_args)
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse the arguments
parsed_args <- parse_args(args)

# Function to check and assign arguments
check_and_assign <- function(arg_name, parsed_args) {
  if (!(arg_name %in% names(parsed_args))) {
    stop(paste("Please provide the", arg_name, "value using the -", arg_name, "argument."))
  }
  return(parsed_args[[arg_name]])
}

# Check and assign dnascreen_run, vcfs_and_bams_dir, sim_dir and bed_file
dnascreen_run <- check_and_assign("dnascreen_run", parsed_args)
vcfs_and_bams_dir <- check_and_assign("vcfs_and_bams_dir", parsed_args)
sim_dir <- check_and_assign("sim_dir", parsed_args)
bed_file <- check_and_assign("bed_file", parsed_args)

# Use regular expressions to extract the desired portion
bed_name <- sub(".*/([^/]+)\\.bed$", "\\1", bed_file)

ref_fa = paste0(sim_dir, "/Homo_sapiens_assembly38.fasta")

# Set the directory path
directory_bam <- paste0(vcfs_and_bams_dir, "/bams_run", dnascreen_run)
directory_vcf <- paste0(vcfs_and_bams_dir, "/vcfs_run", dnascreen_run)

# List all files in the directory with the extension '*.bam'
bam_files <- list.files(directory_bam, pattern = ".bam$", full.names = TRUE)
bai_files <- list.files(directory_bam, pattern = ".bai$", full.names = TRUE)
vcf_files <- list.files(directory_vcf, pattern = ".vcf.gz$", full.names = TRUE)
vcf_filenames <- basename(vcf_files)

extracted <- sub("\\.sorted\\.vcf\\.gz$", "", vcf_filenames)

# Extract sample IDs from bam_files using regex
sample_ids <- sub(".*/([^/]+)\\.hq\\.sorted\\.marked\\.bam", "\\1", bam_files)
sample_bai_ids <- sub(".*/([^/]+)\\.hq\\.sorted\\.marked\\.bam.bai", "\\1", bai_files)

# Check if sample IDs exist in 'extracted' and 'sample_bai_ids'
existing_samples <- sample_ids %in% extracted & sample_ids %in% sample_bai_ids

# Filter bam_files based on existing sample IDs
filtered_bam_files <- bam_files[existing_samples]

filtered_bam_filenames <- basename(filtered_bam_files)
filtered_bam_filenames <- gsub("-", ".", filtered_bam_filenames)

# Ensure each sample is included at least once in 10 sets of 50 samples
set.seed(123)  # For reproducibility
total_samples <- length(filtered_bam_files)
num_sets <- 10
samples_per_set <- 50

# Function to create sets ensuring each sample is included at least once
create_reference_sets <- function(samples, num_sets, samples_per_set) {
  # Calculate number of chunks and chunk sizes
  chunk_sizes <- rep(floor(total_samples / num_sets), num_sets)
  remainder <- total_samples %% num_sets
  chunk_sizes[1:remainder] <- chunk_sizes[1:remainder] + 1
  
  # Create chunks
  chunks <- split(samples, rep(1:num_sets, times = chunk_sizes))
  
  reference_sets <- list()
  for (i in 1:num_sets) {
    set_samples <- unlist(lapply(chunks, function(chunk) sample(chunk, min(samples_per_set / num_sets, length(chunk)), replace = FALSE)))
    set_samples <- sample(set_samples, samples_per_set, replace = FALSE)  # Ensure 50 samples per set
    names(set_samples) <- NULL
    reference_sets[[i]] <- set_samples
  }
  
  return(reference_sets)
}

Rdata_dir <- paste0(sim_dir, "/RData")
reference_sets_file <- paste0(Rdata_dir, "/reference_sets_run", dnascreen_run, ".RData")

if (file.exists(reference_sets_file)) {
  # Load the file if it exists
  load(reference_sets_file)
} else {
  # create reference set if it doesn't exist
  reference_sets <- create_reference_sets(filtered_bam_filenames, 
                                          num_sets, samples_per_set)
  # Save the chosen reference sets to a Rdata file
  save(reference_sets, file = reference_sets_file)
}

# Loading original BAM files
my.count_file <- paste0(Rdata_dir, "/run", dnascreen_run, "_bams_", bed_name, ".RData")
if (file.exists(my.count_file)) {
  # Load the file if it exists
  load(my.count_file)
} else {
  # Run the code to generate and save the file if it does not exist
  my.count <- getBamCounts(bed.file = bed_file,
                           bam.files = filtered_bam_files,
                           include.chr = FALSE,
                           referenceFasta = ref_fa)
  save(my.count, file = my.count_file)
}

my.count.dafr <- as(my.count, 'data.frame')

# Store the results of all.panel@CNV.calls for each file in each run
all_runs_results_list_dnascreen <- list()
ref_sample_lst <- list()

# Run the analysis 10 times with a different reference set
for (run in 1:10) {
  cat("Running analysis iteration:", run, "\n")
  
  # Store the results for the current run
  results_list <- list()
  
  # Random sampling of 50 samples from the filtered_bam_filenames list, excluding the test sample
  dnascreen_samples <- colnames(my.count.dafr[6:length(my.count.dafr)])
  my.ref.samples_for_run <- sample(dnascreen_samples, 50, replace = FALSE)
  ref_sample_lst[[run]] <- my.ref.samples_for_run
  
  # Iterate over each bam file
  for (i in seq(6, ncol(my.count.dafr))) {
    test_filename <- colnames(my.count.dafr)[i]
    
    my.test <- my.count.dafr[[test_filename]]
    
    
    if (sum(my.test == 0) > 10) {
      cat("It looks like the test sample", test_filename, "has more than",
          sum(my.test == 0), "zero reads.",
          "The coverage is too small to perform any meaningful inference so no likelihood will be computed.\n")
      next  # Skip further processing for this test sample
    }
    
    my.ref.samples <- my.ref.samples_for_run[!my.ref.samples_for_run %in% test_filename]
    my.reference.set <- as.matrix(my.count.dafr[, my.ref.samples])
    
    my.choice <- select.reference.set(
      test.counts = my.test,
      reference.counts = my.reference.set,
      bin.length = (my.count.dafr$end - my.count.dafr$start) / 1000,
      n.bins.reduced = 10000
    )
    # Remove duplicates and additional suffix
    my.choice$reference.choice <- unique(gsub("\\.1$", "", my.choice$reference.choice))
    
    my.matrix <- as.matrix(my.count.dafr[, my.choice$reference.choice, drop = FALSE])
    my.reference.selected <- apply(X = my.matrix, MAR = 1, FUN = sum)
    
    all.panel <- new('ExomeDepth',
                     test = my.test,
                     reference = my.reference.selected,
                     formula = 'cbind(test, reference) ~ 1')
    
    all.panel <- CallCNVs(x = all.panel,
                          transition.probability = 10^-4,
                          chromosome = my.count.dafr$chromosome,
                          start = my.count.dafr$start,
                          end = my.count.dafr$end,
                          name = my.count.dafr$exon)
    
    # Check if all.panel@CNV.calls is not empty before adding sample and correlation columns
    if (nrow(all.panel@CNV.calls) > 0) {
      all.panel@CNV.calls$sample <- test_filename
      all.panel@CNV.calls$correlation <- my.choice$summary.stats$correlations[1]
      results_list[[test_filename]] <- all.panel@CNV.calls
    } else {
      # Create a data frame with sample and correlation values if there are no CNV calls
      no_calls_df <- data.frame(
        sample = test_filename,
        correlation = my.choice$summary.stats$correlations[1],
        chromosome = NA,
        start = NA,
        end = NA,
        name = NA,
        type = NA,
        length = NA,
        reads.expected = NA,
        reads.observed = NA,
        reads.ratio = NA,
        stringsAsFactors = FALSE
      )
      results_list[[test_filename]] <- no_calls_df
    }
  }
  
  # Create an empty list to store non-empty data frames
  non_empty_results <- list()
  
  # Iterate over each element in results_list
  for (filename in names(results_list)) {
    # Get the current result
    result <- results_list[[filename]]
    
    # Check if there is a call in the result
    if (!is.data.frame(result) || !all(is.na(result$chromosome))) {
      # Store the non-empty result along with the filename in non_empty_results list
      result$Test_Filename <- filename  # Add filename as a column
      non_empty_results[[filename]] <- result
    }
  }
  
  # Store the non-empty results of the current run in the all_runs_results_list_dnascreen
  all_runs_results_list_dnascreen[[paste0("Run_", run)]] <- non_empty_results
}


# Create an empty list to store the frequency of each unique CNV call across all runs
sample_cnv_calls <- list()

# Iterate over each run in the all_runs_results_list_dnascreen
for (run_name in names(all_runs_results_list_dnascreen)) {
  run_results <- all_runs_results_list_dnascreen[[run_name]]
  
  # Iterate over each sample in the current run's results
  for (sample_name in names(run_results)) {
    cnv_calls <- run_results[[sample_name]]
    
    # Create a unique key based on sample, chr, start, and end for each call
    call_keys <- paste(sample_name, cnv_calls$chromosome, cnv_calls$start, cnv_calls$end, sep = "_")
    
    # Iterate over each call in the current sample
    for (i in 1:nrow(cnv_calls)) {
      call_key <- call_keys[i]
      
      # Check if this call key already exists in the sample_cnv_calls list
      if (call_key %in% names(sample_cnv_calls)) {
        # Increment the frequency count for the call
        sample_cnv_calls[[call_key]]$Frequency <- sample_cnv_calls[[call_key]]$Frequency + 1
        
        # Update the mean values
        sample_cnv_calls[[call_key]]$Mean_Reads_Expected <- 
          round((sample_cnv_calls[[call_key]]$Mean_Reads_Expected * 
                   (sample_cnv_calls[[call_key]]$Frequency - 1) + cnv_calls$reads.expected[i]) / sample_cnv_calls[[call_key]]$Frequency)
        sample_cnv_calls[[call_key]]$Mean_Reads_Observed <- 
          round((sample_cnv_calls[[call_key]]$Mean_Reads_Observed * 
                   (sample_cnv_calls[[call_key]]$Frequency - 1) + cnv_calls$reads.observed[i]) / sample_cnv_calls[[call_key]]$Frequency)
        sample_cnv_calls[[call_key]]$Mean_Reads_Ratio <- 
          round((sample_cnv_calls[[call_key]]$Mean_Reads_Ratio * 
                   (sample_cnv_calls[[call_key]]$Frequency - 1) + cnv_calls$reads.ratio[i]) / sample_cnv_calls[[call_key]]$Frequency, 2)
        sample_cnv_calls[[call_key]]$Mean_BF <- 
          round((sample_cnv_calls[[call_key]]$Mean_BF * 
                   (sample_cnv_calls[[call_key]]$Frequency - 1) + cnv_calls$BF[i]) / sample_cnv_calls[[call_key]]$Frequency, 2)
        
      } else {
        # Initialize the call frequency and store the sample name and CNV call details
        sample_cnv_calls[[call_key]] <- list(
          Frequency = 1,
          Sample = sample_name,
          Chromosome = cnv_calls$chromosome[i],
          Start = cnv_calls$start[i],
          End = cnv_calls$end[i],
          Type = cnv_calls$type[i],
          Nexons = cnv_calls$nexons[i],
          Mean_Reads_Expected = round(cnv_calls$reads.expected[i]),
          Mean_Reads_Observed = round(cnv_calls$reads.observed[i]),
          Mean_Reads_Ratio = round(cnv_calls$reads.ratio[i], 2),
          Mean_BF = round(cnv_calls$BF[i], 2)
        )
      }
    }
  }
}

# Convert the sample_cnv_calls list to a data frame for easier viewing
sample_cnv_calls_df_dnascreen <- do.call(rbind, lapply(sample_cnv_calls, as.data.frame))
sample_cnv_calls_df_dnascreen <- sample_cnv_calls_df_dnascreen[order(sample_cnv_calls_df_dnascreen$Frequency, decreasing = TRUE), ]
sample_cnv_calls_df_dnascreen <- data.frame(sample_cnv_calls_df_dnascreen)

# Create and assign dnascreen run name to 'sample_cnv_calls_df_dnascreen'
# and 'all_runs_results_list_dnascreen'
sample_cnv_calls_var_name <- paste0("sample_cnv_calls_df_dnascreen_run", dnascreen_run)
assign(sample_cnv_calls_var_name, sample_cnv_calls_df_dnascreen)
all_runs_results_list_var_name <- paste0("all_runs_results_list_dnascreen_run", dnascreen_run)
assign(all_runs_results_list_var_name, all_runs_results_list_dnascreen)

# # View the sample CNV calls with their frequencies
# View(sample_cnv_calls_df_dnascreen)

sample_lst <- gsub("\\.hq\\.sorted\\.marked\\.bam$", "", sample_cnv_calls_df_dnascreen$Sample)
sample_lst <- gsub("\\.", "-", sample_lst)
print(sample_lst)

# Write the sample_lst into a text file
writeLines(sample_lst, paste0(sim_dir, "/", "/run", dnascreen_run, "/_samples_with_cnvcalls.txt"))

# saving the sample_cnv_calls_df
write.csv(sample_cnv_calls_df_dnascreen, file = paste0(sim_dir, '/exome_depth_runs_different_ref_samples_results/exome_depth_run', dnascreen_run,'_different_ref_samples.csv'), row.names = FALSE)
