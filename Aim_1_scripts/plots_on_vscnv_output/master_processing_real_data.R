library(ggplot2)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(stringr)
library(GenomicRanges)

plp_cnvs_dir <- "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/clinvar_CNV_list/clinvar_plp"

# List all CSV files ending with "_with_exon.csv"
csv_files <- list.files(path = plp_cnvs_dir, pattern = "*_with_exon.csv", full.names = TRUE)

# Read and combine all CSV files into one data frame
plp_cnvs <- lapply(csv_files, read_and_convert) %>% bind_rows()

# Directories
varseq_dir <- '/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/VarSeq_CNV'
dir_sample_cor <- '/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/sample_correlations'

# Initialize empty data frames
sample_cnv_calls_combined <- data.frame()

# Define dnascreen runs
# dnascreen_runs_real <- c(6, 18, 13, 14, 29, 33)
dnascreen_runs_real <- c(4:35)
# dnascreen_runs_real <- c("_validation_only")

# Process real data and mark outliers
for (dnascreen_run in dnascreen_runs_real) {
  outlier_file_path <- paste0(dir_sample_cor, '/outlier_samples/outlier_IDs_run', dnascreen_run, '.rds')
  # outlier_file_path <- NULL
  
  # Process and mark single-exon data
  sample_cnv_calls_combined <- bind_rows(sample_cnv_calls_combined, process_data_real(dnascreen_run, varseq_dir, outlier_file_path, plp_cnvs, gene_data))
}


# Adjust p-values for multiple testing
adjust_p_values <- function(df, method = "BH") {
  # Adjust p.values using p.adjust
  df <- df %>%
    group_by(run) %>%
    mutate(
      adjusted_p_value = p.adjust(p.value.of.CNV, method = method),
      log_adjusted_p_value = -log(adjusted_p_value, 10)
    )
  return(df)
}

# Apply the multiple-test correction on the combined real data
sample_cnv_calls_combined <- adjust_p_values(sample_cnv_calls_combined, method = "BH")

# Add Number_of_Samples directly to sample_cnv_calls_combined
sample_cnv_calls_combined <- sample_cnv_calls_combined %>%
  group_by(coordinates, Estimated.CN.of.CNV) %>%   # Group by unique CNV calls
  mutate(
    Number_of_Samples = n()                          # Count the number of samples per CNV
  ) %>%
  ungroup() %>%
  arrange(desc(Number_of_Samples)) 

# Ensure Estimated.CN.of.CNV includes all CN values as factors
sample_cnv_calls_combined <- sample_cnv_calls_combined %>%
  mutate(cn = Estimated.CN.of.CNV)

sample_cnv_calls_combined <- sample_cnv_calls_combined %>%
  mutate(data_type = 'Real')

# Below is adding a column that indicates whether sample is sent for acgh validation based on this file below that i sent to Dan.
sent_for_acgh_cnvs <- read.csv("/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/Validations/files_sent_to_dan/cnvs_passed_three_cutoffs.csv")
sent_for_acgh_cnvs$cn <- factor(sent_for_acgh_cnvs$cn, levels = c(0, 1, 2, 3, 4))

# Add validation column based on matching conditions
sample_cnv_calls_combined_check_validation <- sample_cnv_calls_combined %>%
  left_join(sent_for_acgh_cnvs, by = c("sample" = "sample_full_id", "cn", "coordinates")) %>%
  mutate(sent_for_acgh = ifelse(is.na(exon), FALSE, TRUE)) %>%
  select(-starts_with("extra_columns_from_acgh"))  # Remove unnecessary columns from sent_for_acgh_cnvs if needed

sample_cnv_calls_combined$sent_for_acgh <- sample_cnv_calls_combined_check_validation$sent_for_acgh
# Since DNS001303 is requested for aCGH, set it TRUE for sent_for_acgh
sample_cnv_calls_combined[which(sample_cnv_calls_combined$sample_id_trimmed == "DNS001303"), "sent_for_acgh"] <- TRUE

# write.csv(sample_cnv_calls_combined, "varseq_sample_cnv_calls_combined_all_runs_v2.csv")

