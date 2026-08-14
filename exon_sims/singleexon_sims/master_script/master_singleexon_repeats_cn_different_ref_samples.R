library(ExomeDepth)
library(tools)
library(dplyr)

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
directory_bam_test <- paste0(sim_dir, "/singleexon_sims/run", dnascreen_run)
subdirectories_bam_test <- paste0("output_cn", 0:4)

# List all files in the directory with the extension '*.bam'
bam_files_test <- character()
for (subdir in subdirectories_bam_test) {
  dir_path <- file.path(directory_bam_test, subdir)
  files <- list.files(dir_path, pattern = ".bam$", full.names = TRUE)
  bam_files_test <- c(bam_files_test, files)
}

RData_dir <- paste0(sim_dir, "/RData")

# Loading original BAM files
my.count_file <- paste0(RData_dir, "/run", dnascreen_run, "_bams_", bed_name, ".RData")
if (file.exists(my.count_file)) {
  # Load the file if it exists
  load(my.count_file)
} else {
  stop(paste("The file", my.count_file, "does not exist. Stopping script execution."))
}

# Loading 100 modified BAM files
test.count_file <- paste0(RData_dir, "/run", dnascreen_run, "_singleexon_sim_cn_bams_", bed_name, ".RData")
if (file.exists(test.count_file)) {
  # Load the file if it exists
  load(test.count_file)
} else {
  # Run the code to generate and save the file if it does not exist
  test.count <- getBamCounts(bed.file = bed_file,
                             bam.files= bam_files_test,
                             include.chr = FALSE,
                             referenceFasta = ref_fa)
  save(test.count, file = paste0(RData_dir, "/run", dnascreen_run, "_singleexon_sim_cn_bams_", bed_name, ".RData"))
}


my.count.dafr <- as(my.count, 'data.frame')

# Load chosen reference sets
reference_sets_file <- paste0(RData_dir, "/reference_sets_run", dnascreen_run, ".RData")
if (file.exists(reference_sets_file)) {
  # Load the file if it exists
  load(reference_sets_file)
} else {
  stop(paste("The file", reference_sets_file, "does not exist. Stopping script execution."))
}

# Loop 1: Call CNVs for each sample in each cn and store the results

# Store the results of all.panel@CNV.calls for each file in each run
all_runs_results_list_singleexon <- list()
sample_with_no_cnv_calls <- list()
ref_sample_lst <- list()

# Run the analysis 10 times with a different reference set
for (run in 1:10) {
  cat("Running analysis iteration:", run, "\n")
  
  # Store the results for the current run
  results_list <- list()
  
  # Use the pre-defined reference set
  my.ref.samples <- reference_sets[[run]]
  ref_sample_lst[[run]] <- my.ref.samples
  
  # run for each cn
  for (cn in 0:4) {
    # Set the directory paths
    cn_dir <- sprintf('%s/singleexon_sims/run%s/output_cn%d', sim_dir, dnascreen_run, cn)
    selected_regions_dir <- paste0(sim_dir, '/singleexon_sims/run', dnascreen_run,'/selected_regions')
    
    # List all files in the directory with the extension '*.bam'
    cn_files <- list.files(cn_dir, pattern = "\\.bam$", full.names = TRUE)
    
    selected_regions <- read.table(paste0(selected_regions_dir, '/selected_regions_cn', cn, '_sample'), sep = "\t", row.names = NULL)
    colnames(selected_regions) <- c('sample', 'chr', 'start', 'end')
    
    #loops through each sample in selected_regions table
    for (i in 1:nrow(selected_regions)) {
      sel_reg <- selected_regions[i, ]  # Get the current sel_reg as a vector
      sample <- sel_reg['sample'][[1]]
      print(sample)
      test_file <- paste(cn_dir, sprintf("cn%d_%s.bam", cn, sample), sep = "/")
      test_filename <- basename(test_file)
      test_filename <- gsub("-", ".", test_filename)
      my.test <- test.count[[test_filename]]
      
      # my.test <- getBamCounts(bed.file = bed_file,
      #                         bam.files = test_file,
      #                         include.chr = FALSE,
      #                         referenceFasta = ref_fa)
      
      my.ref.samples <- my.ref.samples[!my.ref.samples %in% test_filename]
      my.reference.set <- as.matrix(my.count.dafr[, my.ref.samples])
      
      my.choice <- select.reference.set(test.counts = my.test,
                                        reference.counts = my.reference.set,
                                        bin.length = (my.count.dafr$end - my.count.dafr$start)/1000,
                                        n.bins.reduced = 10000)
      
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
      
      # Remove '.hq.sorted.marked.bam' suffix from each filename in my.choice$reference.choice
      # ref_choice_str <- paste(gsub(".hq.sorted.marked.bam$", "", my.choice$reference.choice), collapse = ",")
      
      # Store the results of all.panel@CNV.calls for the current file
      if (length(all.panel@CNV.calls$chromosome) != 0) {
        results_list[[test_filename]] <- cbind(data.frame(sample = sample, 
                                                          sim_cn = cn, 
                                                          sim_chr = sel_reg$chr, 
                                                          sim_start = sel_reg$start,
                                                          sim_end = sel_reg$end, 
                                                          correlation = my.choice$summary.stats$correlations[1]),
                                               all.panel@CNV.calls)
      }
      else {
        if (cn != 2) {
          empty_run <- data.frame(sample = sample,
                                  run = run,
                                  sim_cn = cn, 
                                  sim_chr = sel_reg$chr, 
                                  sim_start = sel_reg$start,
                                  sim_end = sel_reg$end,
                                  correlation = my.choice$summary.stats$correlations[1])
          
          if (is.null(sample_with_no_cnv_calls[[test_filename]])) {
            sample_with_no_cnv_calls[[test_filename]] <- empty_run
          } else {
            sample_with_no_cnv_calls[[test_filename]] <- rbind(sample_with_no_cnv_calls[[test_filename]], empty_run)
          }
        }
        
      }
    }
  }
  
  # Create an empty list to store non-empty data frames
  non_empty_results <- list()
  
  # Iterate over each element in results_list
  for (filename in names(results_list)) {
    # Get the current result
    result <- results_list[[filename]]
    
    # Check if the result is not an empty data frame
    if (!is.data.frame(result) || (nrow(result) > 0 && ncol(result) > 0)) {
      # Store the non-empty result along with the filename in non_empty_results list
      result$Test_Filename <- filename  # Add filename as a column
      non_empty_results[[filename]] <- result
    }
  }
  # Store the non-empty results of the current run in the all_runs_results_list_singleexon
  all_runs_results_list_singleexon[[paste0("Run_", run)]] <- non_empty_results
}


# explain these dataframes!
sample_cnv_calls <- list()
inaccurate_coordinate_cnvs <- data.frame()
false_positives <- data.frame()


# Loop 2: Check each CNV call stored whether the coordinates match 
# with selected regions
for (i in 1:length(all_runs_results_list_singleexon)) {
  run <- all_runs_results_list_singleexon[[i]]
  for (j in 1:length(run)) {
    cnv_calls <- run[j][[1]]
    
    # check if cnv is within sim region for cn=0,1,3,4
    # check for each cnv call in cnv_calls (even if there is only one, 
    # still important to make this a loop)
    for (k in seq(1, length(cnv_calls$chromosome))) {
      cnv_call <- cnv_calls[k, ]
      sample <- cnv_call$sample
      cn <- cnv_call$sim_cn
      
      if (cnv_call$chromosome == cnv_call$sim_chr) {
        # Calculate the overlap
        overlap_start <- max(cnv_call$start, cnv_call$sim_start)
        overlap_end <- min(cnv_call$end, cnv_call$sim_end)
        overlap_length <- overlap_end - overlap_start
        
        # Calculate the length of the CNV call
        cnv_call_length <- cnv_call$end - cnv_call$start
        
        # Calculate percentage within window
        ED_percentage_within_window <- overlap_length / cnv_call_length
        if (ED_percentage_within_window < 0) {
          ED_percentage_within_window <- 0
        }
      }
      else {
        ED_percentage_within_window <- 0
      }
      
      if (cn != 2) {
        # Create a unique key based on sample, chr, start, and end for each call
        call_key <- paste(sample, cnv_call$chromosome, cnv_call$start, 
                          cnv_call$end, sep = "_")
        # check if call is already in sample_cnv_calls
        if (call_key %in% names(sample_cnv_calls)) {
          # Increment the frequency count for the call
          sample_cnv_calls[[call_key]]$Frequency <- 
            sample_cnv_calls[[call_key]]$Frequency + 1
          
          # Update the mean values and also round ratio and BF to 2 d.p.
          sample_cnv_calls[[call_key]]$Mean_Reads_Expected <- 
            round((sample_cnv_calls[[call_key]]$Mean_Reads_Expected * 
                     (sample_cnv_calls[[call_key]]$Frequency - 1) + cnv_call$reads.expected) / sample_cnv_calls[[call_key]]$Frequency)
          sample_cnv_calls[[call_key]]$Mean_Reads_Observed <- 
            round((sample_cnv_calls[[call_key]]$Mean_Reads_Observed * 
                     (sample_cnv_calls[[call_key]]$Frequency - 1) + cnv_call$reads.observed) / sample_cnv_calls[[call_key]]$Frequency)
          sample_cnv_calls[[call_key]]$Mean_Reads_Ratio <- 
            round((sample_cnv_calls[[call_key]]$Mean_Reads_Ratio * 
                     (sample_cnv_calls[[call_key]]$Frequency - 1) + cnv_call$reads.ratio) / sample_cnv_calls[[call_key]]$Frequency, 2)
          sample_cnv_calls[[call_key]]$Mean_BF <- 
            round((sample_cnv_calls[[call_key]]$Mean_BF * 
                     (sample_cnv_calls[[call_key]]$Frequency - 1) + cnv_call$BF) / sample_cnv_calls[[call_key]]$Frequency, 2)
          sample_cnv_calls[[call_key]]$Mean_Correlation <- 
            round((sample_cnv_calls[[call_key]]$Mean_Correlation * 
                     (sample_cnv_calls[[call_key]]$Frequency - 1) + cnv_call$correlation) / sample_cnv_calls[[call_key]]$Frequency, 2)
        }
        else {
          # fix false positives record for singleexon
          if ((cnv_call$chromosome == cnv_call$sim_chr) &
              (ED_percentage_within_window > 0)) {
            sample_cnv_calls[[call_key]] <- data.frame(sample = sample,
                                                       call_key = call_key,
                                                       Frequency = 1,
                                                       sim_chr = cnv_call$sim_chr,
                                                       sim_start = cnv_call$sim_start,
                                                       sim_end = cnv_call$sim_end,
                                                       sim_length = 
                                                         cnv_call$sim_end - cnv_call$sim_start,
                                                       nexons = cnv_call$nexons,
                                                       ED_length = cnv_call_length,
                                                       ED_chr = cnv_call$chromosome,
                                                       ED_start = cnv_call$start,
                                                       ED_end = cnv_call$end,
                                                       ED_percentage_within_window =
                                                         ED_percentage_within_window,
                                                       Mean_Reads_Expected = 
                                                         cnv_call$reads.expected,
                                                       Mean_Reads_Observed = 
                                                         cnv_call$reads.observed,
                                                       Mean_Reads_Ratio = cnv_call$reads.ratio,
                                                       Mean_BF = cnv_call$BF,
                                                       Mean_Correlation = cnv_call$correlation,
                                                       cn_ED = cn)
          }
          else {
            # records any false positive calls
            false_positives <- rbind(false_positives,
                                     data.frame(sample = sample,
                                                sim_chr = cnv_call$sim_chr,
                                                sim_start = cnv_call$sim_start,
                                                sim_end = cnv_call$sim_end,
                                                sim_length = 
                                                  cnv_call$sim_end - cnv_call$sim_start,
                                                nexons = cnv_call$nexons,
                                                ED_length = cnv_call_length,
                                                ED_chr = cnv_call$chromosome,
                                                ED_start = cnv_call$start,
                                                ED_end = cnv_call$end,
                                                ED_percentage_within_window =
                                                  ED_percentage_within_window,
                                                reads_ratio = cnv_call$reads.ratio,
                                                BF = cnv_call$BF,
                                                cn_ED = cn))
          }
          
          if ((cnv_call$chromosome == cnv_call$sim_chr) &
              (ED_percentage_within_window < 1.0) &
              (ED_percentage_within_window >= 0.9)) {
            inaccurate_coordinate_cnvs <- rbind(inaccurate_coordinate_cnvs, 
                                                data.frame(sample = sample,
                                                           sim_chr = cnv_call$sim_chr,
                                                           sim_start = cnv_call$sim_start,
                                                           sim_end = cnv_call$sim_end,
                                                           sim_length = 
                                                             cnv_call$sim_end - cnv_call$sim_start,
                                                           nexons = cnv_call$nexons,
                                                           ED_length = cnv_call_length,
                                                           ED_chr = cnv_call$chromosome,
                                                           ED_start = cnv_call$start,
                                                           ED_end = cnv_call$end,
                                                           ED_percentage_within_window =
                                                             ED_percentage_within_window,
                                                           reads_ratio = cnv_call$reads.ratio,
                                                           BF = cnv_call$BF,
                                                           cn_ED = cn))
          }
        }
      }
      # if cn = 2, any cnv is a false positive
      else {
        false_positives <- rbind(false_positives,
                                 data.frame(sample = sample,
                                            sim_chr = cnv_call$sim_chr,
                                            sim_start = cnv_call$sim_start,
                                            sim_end = cnv_call$sim_end,
                                            sim_length = 
                                              cnv_call$sim_end - cnv_call$sim_start,
                                            nexons = cnv_call$nexons,
                                            ED_length = cnv_call_length,
                                            ED_chr = cnv_call$chromosome,
                                            ED_start = cnv_call$start,
                                            ED_end = cnv_call$end,
                                            ED_percentage_within_window =
                                              ED_percentage_within_window,
                                            reads_ratio = cnv_call$reads.ratio,
                                            BF = cnv_call$BF,
                                            cn_ED = cn))
      }
    }
  }
}


# Convert the sample_cnv_calls list to a data frame for easier viewing
sample_cnv_calls_df_singleexon <- do.call(rbind, sample_cnv_calls)
sample_cnv_calls_df_singleexon <- sample_cnv_calls_df_singleexon[order(sample_cnv_calls_df_singleexon$Frequency, decreasing = TRUE), ]
sample_cnv_calls_df_singleexon <- data.frame(sample_cnv_calls_df_singleexon)

# Create and assign dnascreen run name to 'sample_cnv_calls_df_singleexon'
# and 'all_runs_results_list_singleexon'
sample_cnv_calls_var_name <- paste0("sample_cnv_calls_df_singleexon_run", dnascreen_run)
assign(sample_cnv_calls_var_name, sample_cnv_calls_df_singleexon)
all_runs_results_list_var_name <- paste0("all_runs_results_list_singleexon_run", dnascreen_run)
assign(all_runs_results_list_var_name, all_runs_results_list_singleexon)

# View(false_positives)
# View(inaccurate_coordinate_cnvs)
# View(sample_cnv_calls_df_singleexon)

# Recording the missed_cnvs
# Convert the sample_with_no_cnv_calls list into a single dataframe
sample_with_no_cnv_calls_df <- do.call(rbind, sample_with_no_cnv_calls)

# Check if sample_with_no_cnv_calls_df is NULL
if (is.null(sample_with_no_cnv_calls_df)) {
  # Create an empty data frame with the same structure
  missed_cnvs <- data.frame(sample = character(0), number_of_missed_runs = integer(0))
} else {
  # Collapse rows by 'sample' and add 'number_of_missed_runs' column
  collapsed_df <- sample_with_no_cnv_calls_df %>%
    group_by(sample) %>%
    summarize(number_of_missed_runs = n()) %>%
    ungroup()
  # Merge the collapsed data with the original data to retain other columns except 'run'
  missed_cnvs <- merge(collapsed_df, sample_with_no_cnv_calls_df, by = "sample")
}

# Merge the collapsed data with the original data to retain other columns except 'run'
missed_cnvs <- merge(collapsed_df, sample_with_no_cnv_calls_df, by = "sample")

# Remove duplicate rows based on 'sample' and select relevant columns
missed_cnvs <- missed_cnvs[!duplicated(missed_cnvs$sample), ]
missed_cnvs <- missed_cnvs[, !(names(missed_cnvs) %in% c("run"))]

write.csv(missed_cnvs, file = paste0(sim_dir, '/singleexon_sims/run', 
                                     dnascreen_run, 
                                     '/plot_generation/missed_cnvs_different_ref_samples.csv'), row.names = FALSE)

write.csv(sample_cnv_calls_df_singleexon, file = paste0(sim_dir, '/singleexon_sims/run', 
                                                       dnascreen_run, 
                                                       '/plot_generation/exome_depth_stats_different_ref_samples.csv'), row.names = FALSE)
write.csv(false_positives, file = paste0(sim_dir, '/singleexon_sims/run', 
                                         dnascreen_run, 
                                         '/plot_generation/false_positives_different_ref_samples.csv'), row.names = FALSE)
write.csv(inaccurate_coordinate_cnvs, file = paste0(sim_dir, '/singleexon_sims/run', 
                                                    dnascreen_run,
                                                    '/plot_generation/inaccurate_coordinate_cnvs_different_ref_samples.csv'), row.names = FALSE)


