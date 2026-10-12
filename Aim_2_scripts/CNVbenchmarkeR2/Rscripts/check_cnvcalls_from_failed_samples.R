library(dplyr)
library(stringr)

# ---- Load failed samples ----
failed_samples <- read.delim("/path/to/project/Aim_1/CNV_calling/check_consent_rescreening/failed_samples_with_bam_sizes.txt", header = FALSE, sep = "\t")

# Clean up into columns
failed_samples <- failed_samples %>%
  mutate(
    Run = as.numeric(str_extract(V1, "\\d+")),  # extract run number
    Sample = str_extract(V2, "DNS-[A-Z0-9-]+_S\\d+"), # extract sample name
    BamSize = str_remove(V3, "BAM size: ")
  ) %>%
  select(Run, Sample, BamSize)

# ---- Check CNV calls against failed samples ----
check_failed <- function(cnv_df, failed_df) {
  cnv_df %>%
    left_join(failed_df, by = c("run" = "Run", "Sample" = "Sample")) %>%
    mutate(FailedSample = !is.na(BamSize))
}

# Apply to each tool
decon_with_failed    <- check_failed(perrun_decon_calls, failed_samples)
clearcnv_with_failed <- check_failed(perrun_clearcnv_calls, failed_samples)
gatk_with_failed     <- check_failed(perrun_gatk_calls, failed_samples)

# ---- Summaries ----
cat("DeCON failed CNV calls:\n")
print(decon_with_failed %>% filter(FailedSample == TRUE))

cat("ClearCNV failed CNV calls:\n")
print(clearcnv_with_failed %>% filter(FailedSample == TRUE))

cat("GATK failed CNV calls:\n")
print(gatk_with_failed %>% filter(FailedSample == TRUE))
