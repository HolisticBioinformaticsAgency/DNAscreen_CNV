library(ExomeDepth)
library(ggplot2)
library(dplyr)
library(tidyr)


# Define the runs you want to analyze
dnascreen_runs <- 4:35  # Add more runs as needed
data_set <- 'real'
bed_name <- '9genes_25bp.fix.sorted'
Rdata_dir <- '/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/exome_depth/RData'

# Container for per-exon data
all_exon_depths <- list()

# Loop over runs
for (dnascreen_run in dnascreen_runs) {
  
  if (data_set == "real") {
    my.count_file <- paste0(
      Rdata_dir, "/run", dnascreen_run,
      "_bams_", bed_name, ".RData"
    )
  } else {
    my.count_file <- paste0(
      Rdata_dir, "/run", dnascreen_run,
      "_multiexon_sim_cn_bams_", bed_name, ".RData"
    )
  }
  
  if (!file.exists(my.count_file)) {
    warning(paste("Missing:", my.count_file))
    next
  }
  
  load(my.count_file)
  
  if (data_set == "simulated") {
    my.count <- test.count
  }
  
  # Extract raw exon counts
  counts <- my.count[, -(1:5)]
  counts <- apply(counts, 2, function(x) as.numeric(trimws(gsub("\"", "", x))))
  counts[is.na(counts)] <- 0
  
  # Exon metadata
  exon_df <- data.frame(
    exon_id = seq_len(nrow(counts)),
    chr = my.count$chromosome,
    start = my.count$start,
    end = my.count$end,
    exon_length_kb = (my.count$end - my.count$start + 1) / 1000
  )
  
  # Sample names
  sample_ids <- colnames(my.count)[6:ncol(my.count)]
  sample_ids <- sub("\\.hq\\.sorted\\.marked\\.bam$", "", sample_ids)
  
  colnames(counts) <- sample_ids
  
  # Convert to long format
  long_df <- as.data.frame(counts) %>%
    mutate(exon_id = exon_df$exon_id) %>%
    pivot_longer(
      cols = -exon_id,
      names_to = "Sample_ID",
      values_to = "Read_Count"
    ) %>%
    left_join(exon_df, by = "exon_id") %>%
    mutate(
      DNAscreen_Run = dnascreen_run
    )
  
  all_exon_depths[[length(all_exon_depths) + 1]] <- long_df
}

combined_exon_depths <- bind_rows(all_exon_depths)

RD_distribution <- ggplot(combined_exon_depths, aes(x = Read_Count)) +
  geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7) +
  scale_x_log10() +
  theme_minimal() +
  labs(
    title = "Distribution of Per-Exon Read Counts",
    x = "Read count per exon (log10)",
    y = "Number of exon–sample observations"
  )

## number of exons that are below 100 reads
low_depth_summary <- combined_exon_depths %>%
  summarise(
    n_total = n(),
    n_below_100 = sum(Read_Count < 100),
    pct_low = 100 * n_below_100 / n_total
  )

low_depth_summary
