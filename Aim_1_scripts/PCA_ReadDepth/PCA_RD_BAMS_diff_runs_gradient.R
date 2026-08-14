library(ggplot2)
library(dplyr)
library(ExomeDepth)
library(plotly)

# Define the runs you want to analyze
dnascreen_runs <- c(4:35)
bed_name <- '9genes_25bp.fix.sorted'
Rdata_dir <- '/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/exome_depth/RData'

# Initialize an empty list to store count matrices
all_counts_list <- list()
all_sample_names <- c()
all_runs <- c()

# Loop over each run to process the data
for (dnascreen_run in dnascreen_runs) {
  my.count_file <- paste0(Rdata_dir, "/run", dnascreen_run, "_bams_", bed_name, ".RData")
  if (file.exists(my.count_file)) {
    # Load the file if it exists
    load(my.count_file)
  } else {
    stop("Error: The file 'my.count_file' does not exist.")
  }
  
  counts <- as.matrix(my.count[, -c(1:5)])  # Remove first 5 columns (non-count data)
  all_counts_list[[dnascreen_run]] <- counts
  
  sample_names <- colnames(counts)
  all_sample_names <- c(all_sample_names, sample_names)
  
  # Assign run number to each sample
  runs <- rep(dnascreen_run, length(sample_names))
  all_runs <- c(all_runs, runs)
}

# Combine all count matrices into a single matrix
combined_counts <- do.call(cbind, all_counts_list)

# Normalize each column (sample) by its total read count
normalized_counts <- apply(combined_counts, 2, function(col) {
  total_reads <- sum(col)
  if (total_reads > 0) {
    (col / total_reads) * 10^3  # Normalize by total read count if it's greater than 0
  } else {
    col  # If total read count is 0, keep the column unchanged
  }
})

# Perform PCA
pca_results <- prcomp(t(normalized_counts), scale. = TRUE)

# Extract percentage of variance explained by PC1 and PC2
explained_variance <- pca_results$sdev^2 / sum(pca_results$sdev^2)
PC1_var <- round(explained_variance[1] * 100, 2)
PC2_var <- round(explained_variance[2] * 100, 2)

# Create a data frame for plotting
pca_df <- as.data.frame(pca_results$x)
pca_df$Sample <- all_sample_names
pca_df$Run <- factor(all_runs, levels = dnascreen_runs)

# Dynamically create a title by pasting the dnascreen_runs
plot_title <- paste("PCA of Samples Based on Read Counts for Runs 4-35")

# Plot PCA with variance explained in axis labels
ggplot(pca_df, aes(x = PC1, y = PC2, color = as.numeric(Run))) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = c("#d0e1f9", "#74a9cf", "#0570b0", "#034e7b")) +  # Light to dark blue shades
  labs(
    title = plot_title,
    x = paste0("Principal Component 1 (", PC1_var, "% variance)"),
    y = paste0("Principal Component 2 (", PC2_var, "% variance)"),
    color = "Run"
  ) +
  theme_minimal()


# For personal visualisation purposes
# Plot PCA with discrete colors for runs
ggplot_discrete <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Run)) +
  geom_point(size = 3) +
  scale_color_manual(values = scales::hue_pal()(length(dnascreen_runs))) +  # Distinct colors for each run
  labs(
    title = plot_title,
    x = paste0("Principal Component 1 (", PC1_var, "% variance)"),
    y = paste0("Principal Component 2 (", PC2_var, "% variance)"),
    color = "Run"
  ) +
  theme_minimal()

# Convert ggplot to plotly for interactivity
plotly_discrete <- ggplotly(ggplot_discrete)

# Save the interactive plot as an HTML file
htmlwidgets::saveWidget(plotly_discrete, 
                        "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/PCA_CNV_detection/PCA_plot_runs_4-35.html")
