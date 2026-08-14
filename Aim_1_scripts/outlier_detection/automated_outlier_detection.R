# Libraries
library(ExomeDepth)
library(tools)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(dbscan)
library(pracma)
library(plotly)
library(htmlwidgets)

# Define the runs you want to analyze
# Vector of runs for looping
dnascreen_runs <- c(4:35)

# Other variables remain the same
bed_name <- '9genes_25bp.fix.sorted'
Rdata_dir <- '/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/exome_depth/RData'
dir_sample_cor <- '/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/sample_correlations'


# Loop over each run in dnascreen_runs
for (dnascreen_run in dnascreen_runs) {
  
  # Define the file path for the current run
  my.count_file <- paste0(Rdata_dir, "/run", dnascreen_run, "_bams_", bed_name, ".RData")
  
  # Check if the file exists, and load it, or skip if not
  if (file.exists(my.count_file)) {
    load(my.count_file)
  } else {
    warning(paste("Warning: The file for run", dnascreen_run, "does not exist. Skipping this run."))
    next  # Skip to the next run if the file doesn't exist
  }
  
  # Transpose and clean the counts matrix
  counts <- t(my.count)[-c(1:5), ]
  counts <- apply(counts, 2, function(x) as.numeric(trimws(gsub("\"", "", x))))
  counts[is.na(counts)] <- 0
  
  # Calculate feature lengths and normalize counts by RPKM
  feature_lengths <- (my.count$end - my.count$start + 1) / 1000
  normalized_counts <- t(apply(counts, 1, function(row) {
    total_reads <- sum(row)
    if (total_reads > 0) {
      rpkm <- ((row / total_reads) / feature_lengths) * 10^6
      rpkm
    } else {
      row
    }
  }))
  
  # Calculate sample-sample correlation matrix
  cor_matrix <- cor(t(normalized_counts))
  cor_matrix_no_diag <- cor_matrix
  diag(cor_matrix_no_diag) <- NA
  
  # Extract sample IDs and process
  sample_ids <- colnames(my.count[6:length(my.count)])
  sample_ids <- gsub("\\.hq\\.sorted\\.marked\\.bam$", "", sample_ids)
  sample_ids <- gsub("\\.", "-", sample_ids)
  sample_ids <- sub(".*-(DNS[0-9]+)_.*", "\\1", sample_ids)
  
  # Save the correlation matrix without diagonal
  saveRDS(cor_matrix_no_diag, file = paste0(dir_sample_cor, '/sample_cor_matrix/cor_matrix_no_diag_run', dnascreen_run, '.rds'))
  
  # Compute the maximum correlation values
  max_correlation_values <- apply(cor_matrix_no_diag, 2, function(col) max(col, na.rm = TRUE))
  names(max_correlation_values) <- sample_ids
  
  # Save a histogram of correlation values
  histogram_plot <- ggplot(data.frame(max_correlation_values), aes(x = max_correlation_values)) +
    geom_histogram(bins = 10, fill = "blue", color = "black") +
    labs(title = paste0("Histogram of Correlation Values for Run ", dnascreen_run),
         x = "Correlation", y = "Frequency") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"))
  
  ggsave(filename = paste0(dir_sample_cor, '/histograms/histogram_run', dnascreen_run, '.png'),
         plot = histogram_plot, width = 8, height = 6)
  
  # Save the max correlation values as CSV
  max_correlation_values_df <- data.frame(
    sample = names(max_correlation_values),
    max_correlation_value = max_correlation_values
  )
  write.csv(max_correlation_values_df, paste0(dir_sample_cor, '/max_correlation_values/max_correlation_values_run', dnascreen_run, '.csv'))
  
  # Perform MDS and DBSCAN
  dissimilarity <- 1 - cor_matrix
  mds <- cmdscale(dissimilarity, k = 2)
  plot_data <- data.frame(X = mds[, 1], Y = mds[, 2], Sample_ID = sample_ids)
  
  # Define clustering parameters and perform DBSCAN
  minPts <- 5
  distances <- dist(mds)
  k_distances <- apply(as.matrix(distances), 1, function(row) sort(row)[minPts])
  k_distances_sorted <- sort(k_distances)
  
  # Determine the elbow point for epsilon
  second_derivative <- diff(diff(k_distances_sorted))
  elbow_index <- which.max(second_derivative) + 1
  epsilon <- k_distances_sorted[elbow_index]
  
  # Apply DBSCAN
  clustering <- dbscan(mds, eps = epsilon, minPts = minPts)
  plot_data$Cluster <- as.factor(clustering$cluster)
  
  # Identify and save outliers
  outliers <- subset(plot_data, Cluster == 0)
  saveRDS(outliers$Sample_ID, file = paste0(dir_sample_cor, '/outlier_samples/outlier_IDs_run', dnascreen_run, '.rds'))
  
  # Plot MDS with clusters and outliers
  p <- ggplot(plot_data, aes(x = X, y = Y, label = Sample_ID, color = Cluster)) +
    geom_point(size = 3) +
    geom_text_repel(data = outliers, aes(label = Sample_ID), size = 3, segment.size = 0.2, segment.color = "black") +
    labs(title = paste0("Sample-sample correlation distance for DNA Screen Run ", dnascreen_run)) +
    theme_minimal() +
    theme(legend.position = "right") +
    scale_color_discrete(name = "Cluster")
  
  # Convert ggplot object to Plotly and save as HTML
  plotly_p <- ggplotly(p)
  saveRDS(plotly_p, file = paste0(dir_sample_cor, '/MDS_plots/MDS_plot_run', dnascreen_run, '.rds'))
  saveWidget(plotly_p, file = paste0(dir_sample_cor, '/MDS_plots/MDS_plot_run', dnascreen_run, '.html'), selfcontained = TRUE)
}
