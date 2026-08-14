library(ggplot2)
library(dplyr)

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
  
  counts <- as.matrix(my.count[, -c(1:5)])
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

# Create a data frame for plotting
pca_df <- as.data.frame(pca_results$x)
pca_df$Sample <- all_sample_names
pca_df$Run <- factor(all_runs, levels = dnascreen_runs)

# Calculate mean PC1 and PC2 for each run
mean_pc_df <- pca_df %>%
  group_by(Run) %>%
  summarize(Mean_PC1 = mean(PC1), Mean_PC2 = mean(PC2))

# Perform k-means clustering with the chosen number of clusters (e.g., k = 3)
set.seed(123)  # For reproducibility
k <- 3 # Choose the number of clusters based on the elbow plot
kmeans_result <- kmeans(mean_pc_df[, c("Mean_PC1", "Mean_PC2")], centers = k, nstart = 10)

# Convert the cluster to numeric or explicitly define labels (this will avoid the "a" issue)
mean_pc_df$Cluster <- factor(kmeans_result$cluster, labels = paste0("Cluster ", 1:k))

# Find the closest two runs to the centroid of each cluster
closest_runs <- mean_pc_df %>%
  group_by(Cluster) %>%
  mutate(
    Centroid_PC1 = mean(Mean_PC1),
    Centroid_PC2 = mean(Mean_PC2),
    Distance_to_Centroid = euclidean_distance(Mean_PC1, Mean_PC2, Centroid_PC1, Centroid_PC2)
  ) %>%
  arrange(Distance_to_Centroid) %>%
  slice_head(n = 2)

ggplot(mean_pc_df, aes(x = Mean_PC1, y = Mean_PC2, color = Cluster)) +
  geom_point(size = 4) +  # Points colored by cluster
  geom_text(aes(label = Run), vjust = -1, hjust = 0.5, show.legend = FALSE) +  # Add run labels for all points
  theme_minimal() +
  labs(title = "K-means Clustering of Runs 4-35 Based on Mean PC1 and PC2",
       x = "Mean PC1",
       y = "Mean PC2") +
  theme(legend.title = element_blank()) +
  geom_text(data = closest_runs, aes(x = Mean_PC1, y = Mean_PC2, label = Run), 
            vjust = -1, hjust = 0.5, color = 'red', fontface = "bold", size = 4, inherit.aes = FALSE) +  # Closest runs with larger, bold text
  scale_color_manual(values = c("Cluster 1" = "#034e7b", "Cluster 2" = "#74a9cf", "Cluster 3" = "#d0e1f9"))  # Define custom colors

