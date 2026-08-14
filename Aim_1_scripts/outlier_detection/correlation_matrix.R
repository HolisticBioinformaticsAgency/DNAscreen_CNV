# Load necessary libraries
library(ExomeDepth)
library(tools)
library(ggplot2)
library(ggrepel)
library(plotly)
library(crosstalk)
library(dbscan)

# Define the runs you want to analyze
dnascreen_runs <- 4:35  # Add more runs as needed
data_set <- 'real'
bed_name <- '9genes_25bp.fix.sorted'
Rdata_dir <- '/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/exome_depth/RData'

# Initialize an empty list to store data
all_plot_data <- list()

# Iterate through each run
for (dnascreen_run in dnascreen_runs) {
  
  # Define the file path for the current run
  if (data_set == 'real') {
    my.count_file <- paste0(Rdata_dir, "/run", dnascreen_run, "_bams_", bed_name, ".RData")
  } else if (data_set == 'simulated') {
    my.count_file <- paste0(Rdata_dir, "/run", dnascreen_run,
                            "_multiexon_sim_cn_bams_", bed_name, ".RData")
  }
  
  # Skip if file doesn't exist
  if (!file.exists(my.count_file)) {
    warning(paste("Warning: The file for run", dnascreen_run, "does not exist. Skipping this run."))
    next
  }
  
  # Load the file
  load(my.count_file)
  
  # Use appropriate data for real or simulated
  if (data_set == 'simulated') {
    my.count <- test.count
  }
  
  # Transpose the counts matrix
  counts <- t(my.count)[-c(1:5), ]
  
  # Clean and convert counts matrix to numeric
  counts <- apply(counts, 2, function(x) as.numeric(trimws(gsub("\"", "", x))))
  
  # Check for NA values and replace with 0
  counts[is.na(counts)] <- 0
  
  # Calculate feature lengths (assuming BED file contains start and end coordinates)
  feature_lengths <- (my.count$end - my.count$start + 1)/1000
  
  # Normalize counts matrix using RPKM
  normalized_counts <- t(apply(counts, 1, function(row) {
    total_reads <- sum(row)
    if (total_reads > 0) {
      rpkm <- ((row / total_reads)/ feature_lengths)*10^6  # RPKM calculation
      rpkm
    } else {
      row
    }
  }))
  
  # Calculate correlation matrix
  cor_matrix <- cor(t(normalized_counts))
  
  # Convert correlation to dissimilarity (distance)
  dissimilarity <- 1 - cor_matrix
  
  # Apply MDS to reduce dimensionality to 2D
  mds <- cmdscale(dissimilarity, k = 2)
  
  # Extract sample IDs
  sample_id_cols <- colnames(my.count)[6:length(my.count)]
  sample_id_cols <- sub("\\.hq\\.sorted\\.marked\\.bam$", "", sample_id_cols)
  
  # Create data frame for plotting
  plot_data <- data.frame(
    X = mds[, 1],
    Y = mds[, 2],
    Sample_ID = sample_id_cols,
    DNAscreen_Run = as.factor(dnascreen_run)
  )
  
  # Append to the list
  all_plot_data[[length(all_plot_data) + 1]] <- plot_data
}

# Combine all runs' data into a single data frame
combined_plot_data <- do.call(rbind, all_plot_data)

# Shareable dataframe for crosstalk
shared_data <- SharedData$new(combined_plot_data, key = ~DNAscreen_Run)

# Generate the plotly correlation plot
interactive_plot <- plot_ly(
  shared_data,
  x = ~X,
  y = ~Y,
  type = 'scatter',
  mode = 'markers',
  color = ~DNAscreen_Run,
  colors = colorRampPalette(c("lightblue", "darkblue"))(length(dnascreen_runs)),
  text = ~paste("Sample ID:", Sample_ID, "<br>Run:", DNAscreen_Run),
  hoverinfo = "text"
) %>%
  layout(
    title = "MDS Plot of Samples (Filter by Run)",
    xaxis = list(title = "Dimension 1"),
    yaxis = list(title = "Dimension 2")
  ) %>%
  highlight(
    on = "plotly_click",
    off = "plotly_doubleclick",
    dynamic = TRUE,
    color = "red"
  )

# Print the interactive plot
interactive_plot
