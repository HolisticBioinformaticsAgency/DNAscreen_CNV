library(ggplot2)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(stringr)
library(GenomicRanges)

# Directories
clinvar_exon_sim_dir <- '/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/simulation_hom_deletion/experiment_onsimulation_using_clinvar'
varseq_dir <- '/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/VarSeq_CNV'
dir_sample_cor <- '/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/sample_correlations'

# Initialize empty data frames
sim_recalled_clinvar_table_combined <- data.frame()
sim_false_positives_clinvar_table_combined <- data.frame()
sim_missed_clinvar_table_combined <- data.frame()
sim_recalled_after_sample_shuffling_clinvar_table_combined <- data.frame()

# Define dnascreen runs
dnascreen_runs_sim <- c(6, 18, 13, 14, 29, 33)
# dnascreen_runs_sim <- c(6)

# Process simulated data
# Note there are no cn=2 simulations for clinvar
for (dnascreen_run in dnascreen_runs_sim) {
  sim_recalled_clinvar_table_combined <- bind_rows(sim_recalled_clinvar_table_combined, process_data_sim_clinvar(dnascreen_run, clinvar_exon_sim_dir, "recalled"))
}

sim_recalled_clinvar_table_combined <- sim_recalled_clinvar_table_combined %>%
  mutate(data_type = 'Simulated',
         within_exon = ifelse(
           is.na(overlap_with_start_name) | is.na(overlap_with_end_name), 
           FALSE, 
           overlap_with_start_name == overlap_with_end_name))
sim_recalled_clinvar_table_combined$log_adjusted_p_value <- sim_recalled_clinvar_table_combined$log_p_value
sim_recalled_clinvar_table_combined$Avg.Read.Ratio.of.CNV <- sim_recalled_clinvar_table_combined$Read.Ratio.of.Simulated.CNV

# Getting missed CNVs
for (dnascreen_run in dnascreen_runs_sim) {
  sim_missed_clinvar_table_combined <- bind_rows(sim_missed_clinvar_table_combined, process_data_sim_clinvar(dnascreen_run, clinvar_exon_sim_dir, "missed_cnvs"))
}

sim_missed_clinvar_table_combined <- sim_missed_clinvar_table_combined %>%
  mutate(data_type = 'Missed',
         within_exon = ifelse(
           is.na(overlap_with_start_name) | is.na(overlap_with_end_name), 
           FALSE, 
           overlap_with_start_name == overlap_with_end_name))
sim_missed_clinvar_table_combined$log_p_value <- rep(0, nrow(sim_missed_clinvar_table_combined))
sim_missed_clinvar_table_combined$log_adjusted_p_value <- rep(0, nrow(sim_missed_clinvar_table_combined))
sim_missed_clinvar_table_combined$Avg.Read.Ratio.of.CNV <- sim_missed_clinvar_table_combined$Read.Ratio.of.Simulated.CNV


# getting false_positives
for (dnascreen_run in dnascreen_runs_sim) {
  sim_false_positives_clinvar_table_combined <- bind_rows(sim_false_positives_clinvar_table_combined, process_data_sim_clinvar(dnascreen_run, clinvar_exon_sim_dir, "false_positive"))
}
sim_false_positives_clinvar_table_combined <- sim_false_positives_clinvar_table_combined %>%
  mutate(data_type = 'False Positive',
         within_exon = ifelse(
           is.na(overlap_with_start_name) | is.na(overlap_with_end_name), 
           FALSE, 
           overlap_with_start_name == overlap_with_end_name))
sim_false_positives_clinvar_table_combined$log_adjusted_p_value <- sim_false_positives_clinvar_table_combined$log_p_value
sim_false_positives_clinvar_table_combined$Avg.Read.Ratio.of.CNV <- sim_false_positives_clinvar_table_combined$Read.Ratio.of.Simulated.CNV

# getting missed CNVs that were recalled on the second try
# for (dnascreen_run in dnascreen_runs_sim) {
#   sim_recalled_after_sample_shuffling_clinvar_table_combined <- bind_rows(sim_recalled_after_sample_shuffling_clinvar_table_combined, process_data_sim_clinvar(dnascreen_run, clinvar_exon_sim_dir, "recalled_after_sample_shuffling"))
# }
# sim_recalled_after_sample_shuffling_clinvar_table_combined <- sim_recalled_after_sample_shuffling_clinvar_table_combined %>%
#   mutate(data_type = 'Recalled after sample shuffling')
# sim_recalled_after_sample_shuffling_clinvar_table_combined$log_adjusted_p_value <- sim_recalled_after_sample_shuffling_clinvar_table_combined$log_p_value


