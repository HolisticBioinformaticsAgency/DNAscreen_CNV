library(dplyr)
library(readr)
library(stringr)

read_all_hsmetrics <- function(hsdir) {
  files <- list.files(hsdir, pattern = "\\.HSmetrics\\.txt$", full.names = TRUE)
  
  all_metrics <- lapply(files, function(f) {
    skip_lines <- grep("^BAIT_SET", readLines(f)) - 1
    read_tsv(f, skip = skip_lines) %>%
      mutate(Sample = str_remove(basename(f), "\\.HSmetrics\\.txt$"))
  }) %>% bind_rows()
  
  return(all_metrics)
}

# Usage
hsdir <- "/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/hsmetrics"
all_metrics <- read_all_hsmetrics(hsdir)

# Quick summary: % on-target per barcode
summary_stats <- all_metrics %>%
  select(Sample, TOTAL_READS, ON_TARGET_BASES, PF_BASES, PCT_USABLE_BASES_ON_TARGET, PCT_TARGET_BASES_10X, PCT_TARGET_BASES_20X, PCT_TARGET_BASES_30X, PCT_TARGET_BASES_40X, PCT_TARGET_BASES_50X, PCT_TARGET_BASES_100X, AT_DROPOUT, GC_DROPOUT) %>%
  mutate(Pct_on_target = ON_TARGET_BASES / PF_BASES * 100) %>%
  # keep rows where at least one metric column is not NA
  filter(!if_all(c(TOTAL_READS, ON_TARGET_BASES, PF_BASES, PCT_USABLE_BASES_ON_TARGET), is.na))

# View or save
print(summary_stats)
write_csv(summary_stats, file = file.path(hsdir, "hsmetrics_per_barcode_summary.csv"))
