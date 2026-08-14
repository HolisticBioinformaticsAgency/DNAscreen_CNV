library(ggplot2)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(stringr)
library(GenomicRanges)

# Read the gene coordinates
gene_data <- read.table("dnascreen_genes.txt", header = FALSE, col.names = c("chr", "start", "end", "Gene.Names"))
gene_data$chr <- str_replace(gene_data$chr, "^chr", "")

# Function to assign gene names to CNV calls
assign_gene_names <- function(cnv_calls, gene_data) {
  # Create GenomicRanges for CNV calls
  cnv_gr <- GRanges(
    seqnames = cnv_calls$VS_Call_chr,
    ranges = IRanges(
      start = cnv_calls$VS_Call_start,
      end = cnv_calls$VS_Call_end
    )
  )
  
  # Create GenomicRanges for gene data
  gene_gr <- GRanges(
    seqnames = gene_data$chr,
    ranges = IRanges(start = gene_data$start, end = gene_data$end),
    Gene.Names = gene_data$Gene.Names
  )
  
  # Find overlaps between CNV regions and gene regions
  overlaps <- findOverlaps(cnv_gr, gene_gr)
  
  # Initialize the 'Gene.Names' column with NA
  cnv_calls$Gene.Names <- NA
  
  # Assign gene names to the overlapping CNV regions
  cnv_calls$Gene.Names[queryHits(overlaps)] <- gene_gr$Gene.Names[subjectHits(overlaps)]
  
  return(cnv_calls)
}


# Helper function to check if two CNV intervals overlap
check_overlap <- function(chr1, start1, end1, chr2, start2, end2) {
  # Ensure the chromosome is the same and there is an overlap in the start-end range
  chr1 == chr2 && (start1 <= end2 && end1 >= start2)
}


# Function to read CSV and ensure Exon_Info is treated as character
read_and_convert <- function(file) {
  df <- read.csv(file)
  
  # Convert 'Exon_Info' to character if it exists
  if ("Exon_Info" %in% names(df)) {
    df$Exon_Info <- as.character(df$Exon_Info)
  }
  
  return(df)
}


# Function to process simulated single-exon CNVs
process_data_sim_single <- function(dnascreen_run, sim_dir, data_type) {
  # Read in the simulation table
  if (data_type == 'recalled'){
    sim_table <- read.csv(paste0(sim_dir, '/run', 
                                 dnascreen_run, 
                                 '/plot_generation/Recalled_VS_CNV_Calls_metrics.csv'))
    sim_table$recalled <- TRUE
  } 
  if (data_type == 'false_positive'){
    sim_table <- read.csv(paste0(sim_dir, '/run', 
                                 dnascreen_run, 
                                 '/plot_generation/varseq_tables/false_positives_varseq.csv'))
    sim_table$recalled <- FALSE
  }
  if (data_type == 'cn2'){
    sim_table <- read.csv(paste0(sim_dir, '/run', 
                                 dnascreen_run, 
                                 '/plot_generation/cn_2_sims.csv'))
    sim_table$recalled <- FALSE
  }
  if (data_type == 'missed_cnvs'){
    sim_table <- read.csv(paste0(sim_dir, '/run', 
                                 dnascreen_run, 
                                 '/plot_generation/varseq_tables/missed_cnvs_with_varseq_calls.csv'))
    sim_table$recalled <- FALSE
  }
  if (data_type == 'recalled_after_sample_shuffling'){
    sim_table_list <- list()
    suffix_number <- 1
    
    while (TRUE) {
      # Construct the file path based on the current suffix number
      suffix <- ifelse(suffix_number == 1, "_missed", paste0("_missed_", suffix_number))
      file_path <- paste0(sim_dir, "/run", dnascreen_run, suffix, "/plot_generation/Recalled_VS_CNV_Calls_metrics.csv")
      
      # Check if the file exists
      if (!file.exists(file_path)) {
        message("File not found: ", file_path, ". Terminating the loop.")
        break
      }
      
      # Read the CSV file and mark it as recalled
      sim_table <- read.csv(file_path)
      sim_table$recalled <- TRUE
      
      # Add the data to the list (optional, if you need to store it)
      sim_table_list[[suffix]] <- sim_table
      
      # Increment the suffix number
      suffix_number <- suffix_number + 1
    }
    
    # Optionally combine all data frames into one if needed
    if (length(sim_table_list) > 0) {
      sim_table <- do.call(rbind, sim_table_list)
    }
    
  }
  
  # Extract the relevant part of the sample ID (e.g., "DNS001451" from "DNS-XTHS-0018-B09-DNS001451_S258")
  sim_table <- sim_table %>%
    mutate(sample_id_trimmed = str_extract(sample, "DNS\\d+"))  # Extracts the DNS ID pattern
  
  # Log transformation for simulation data and filtering based on exon type
  sim_table <- sim_table %>%
    mutate(
      # Remove 'chr' from VS_Call_chr
      VS_Call_chr = str_replace(VS_Call_chr, "^chr", ""),
      # Create a new 'coordinates' column by combining VS_Call_chr, VS_Call_start, and VS_Call_end
      coordinates = paste0(VS_Call_chr, ":", VS_Call_start, "-", VS_Call_end),
      simulated_depth = Target.Mean.Depth.of.Simulated.CNV,
      simulated_rr = Read.Ratio.of.Simulated.CNV,
      simulated_z_score = Z.Score.of.Simulated.CNV,
      log_p_value = ifelse(is.nan(p.value.of.CNV), NaN, -log(p.value.of.CNV, 10)),  # Log-transform the p-value
      cn = factor(cn, levels = c(0, 1, 2, 3, 4)),
      Estimated.CN.of.CNV = factor(Estimated.CN.of.CNV, levels = c(0, 1, 2, 3, 4)),
      run = dnascreen_run,
      Very.High.Sensitivity = as.logical(Very.High.Sensitivity),
      High.Sensitivity = as.logical(High.Sensitivity),
      Balanced = as.logical(Balanced),
      High.Precision = as.logical(High.Precision),
      Very.High.Precision = as.logical(Very.High.Precision),
      # Create the tooltip text with rounding
      text = paste("Sample:", sample, "<br>",
                   "Coordinates:", coordinates, "<br>",
                   "P-Value:", round(p.value.of.CNV, 3), "<br>",
                   "Called at precision levels:", Precision.Level.of.call., "<br>",
                   "Simulated CN:", cn, "<br>",
                   "Estimated CN:", Estimated.CN.of.CNV, "<br>",
                   "Mean Depth of Simulated CNV:", round(simulated_depth, 2), "<br>",
                   "Z-Score of Simulated CNV:", round(simulated_z_score, 2), "<br>",
                   "Read Ratio of Simulated CNV:", round(simulated_rr, 3), "<br>",
                   "Mean Depth of VS-called CNV:", round(VS.Target.Mean.Depth.of.CNV, 2), "<br>",
                   "Z-Score of VS-called CNV:", round(VS.Z.Score.of.CNV, 2), "<br>",
                   "Read Ratio of VS-called CNV:", round(VS.Read.Ratio.of.CNV, 3), "<br>",
                   "Run:", dnascreen_run, "<br>",
                   "Percent Difference:", round(Percent.Difference, 2), "<br>",
                   "Recalled:", recalled, "<br>")
    )
  
  sim_table
}

process_data_sim_clinvar <- function(dnascreen_run, sim_dir, data_type) {
  # Read in the simulation table
  if (data_type == 'recalled'){
    sim_table <- read.csv(paste0(sim_dir, '/run', 
                                 dnascreen_run, 
                                 '/plot_generation/Recalled_VS_CNV_Calls_metrics.csv'))
    sim_table$recalled <- TRUE
  } 
  if (data_type == 'false_positive'){ 
    sim_table <- read.csv(paste0(sim_dir, '/run', 
                                 dnascreen_run, 
                                 '/plot_generation/varseq_tables/false_positives_varseq.csv'))
    sim_table$recalled <- FALSE
  }
  if (data_type == 'missed_cnvs'){
    sim_table <- read.csv(paste0(sim_dir, '/run', 
                                 dnascreen_run, 
                                 '/plot_generation/varseq_tables/missed_cnvs_with_varseq_calls.csv'))
    sim_table$recalled <- FALSE
  }
  if (data_type == 'recalled_after_sample_shuffling'){
    sim_table <- read.csv(paste0(sim_dir, '/run', 
                                 dnascreen_run, '_missed',
                                 '/plot_generation/Recalled_VS_CNV_Calls_metrics.csv'))
    sim_table$recalled <- TRUE
  }
  
  # Extract the relevant part of the sample ID
  sim_table <- sim_table %>%
    mutate(sample_id_trimmed = str_extract(sample, "DNS\\d+"))
  
  # Log transformation for simulation data and filtering based on exon type
  sim_table <- sim_table %>%
    mutate(
      # Remove 'chr' from VS_Call_chr
      VS_Call_chr = str_replace(VS_Call_chr, "^chr", ""),
      # Create a new 'coordinates' column by combining VS_Call_chr, VS_Call_start, and VS_Call_end
      coordinates = paste0(VS_Call_chr, ":", VS_Call_start, "-", VS_Call_end),
      orig_depth = Target.Mean.Depth.of.Original.CNV.region,
      orig_rr = Read.Ratio.of.Original.CNV.region,
      orig_z_score = Z.Score.of.Original.CNV.region,
      simulated_depth = Target.Mean.Depth.of.Simulated.CNV,
      simulated_rr = Read.Ratio.of.Simulated.CNV,
      simulated_z_score = Z.Score.of.Simulated.CNV,
      log_p_value = ifelse(is.nan(p.value.of.CNV), NaN, -log(p.value.of.CNV, 10)),
      cn = factor(cn, levels = c(0, 1, 2, 3, 4)),
      Estimated.CN.of.CNV = factor(Estimated.CN.of.CNV, levels = c(0, 1, 2, 3, 4)),
      run = dnascreen_run,
      Very.High.Sensitivity = as.logical(Very.High.Sensitivity),
      High.Sensitivity = as.logical(High.Sensitivity),
      Balanced = as.logical(Balanced),
      High.Precision = as.logical(High.Precision),
      Very.High.Precision = as.logical(Very.High.Precision),
      # Create the tooltip text with rounded values
      text = paste("Sample:", sample, "<br>",
                   "Coordinates:", coordinates, "<br>",
                   "P-Value:", round(p.value.of.CNV, 4), "<br>",
                   "Called at precision levels:", Precision.Level.of.call., "<br>",
                   "Simulated CN:", cn, "<br>",
                   "Estimated CN:", Estimated.CN.of.CNV, "<br>",
                   "Number of exons within simulated window", count_exons_within_simulated_window, "<br>",
                   "Mean Depth of Original CNV region:", round(orig_depth, 2), "<br>",
                   "Z-Score of Original CNV region:", round(orig_z_score, 2), "<br>",
                   "Read Ratio of Original CNV region:", round(orig_rr, 2), "<br>",
                   "Mean Depth of Simulated CNV:", round(simulated_depth, 2), "<br>",
                   "Z-Score of Simulated CNV:", round(simulated_z_score, 2), "<br>",
                   "Read Ratio of Simulated CNV:", round(simulated_rr, 2), "<br>",
                   "Mean Depth of VS-called CNV:", round(VS.Target.Mean.Depth.of.CNV, 2), "<br>",
                   "Z-Score of VS-called CNV:", round(VS.Z.Score.of.CNV, 2), "<br>",
                   "Read Ratio of VS-called CNV:", round(VS.Read.Ratio.of.CNV, 2), "<br>",
                   "Run:", dnascreen_run, "<br>",
                   "Percent Difference:", round(Percent.Difference, 2), "<br>",
                   "Recalled:", recalled, "<br>",
                   "Total overlap length with overlapping exons:", round(exon_overlap_length, 2), "<br>",
                   "Percentage total length of overlapping exons that overlaps with simulated window:",
                   round(exon_overlap_percentage, 2), "<br>")
    )
  
  sim_table
}


# Function to process real data and mark outliers (only if specified)
process_data_real <- function(dnascreen_run, varseq_dir, outlier_file_path = NULL, plp_cnvs, gene_data) {
  # Read in the CNV calls
  sample_cnv_calls <- read.csv(paste0(varseq_dir, '/varseq_tables_on_real_data', '/run', dnascreen_run, '/CNV_in_each_sample_run', dnascreen_run, '.csv'))
  
  # Extract the relevant part of the sample ID (e.g., "DNS001451" from "DNS-XTHS-0018-B09-DNS001451_S258")
  sample_cnv_calls <- sample_cnv_calls %>%
    mutate(sample_id_trimmed = str_extract(sample, "DNS\\d+"))  # Extracts the DNS ID pattern
  
  # Load outliers only if outlier_file_path is provided
  if (!is.null(outlier_file_path) && outlier_file_path != "") {
    outliers <- readRDS(outlier_file_path)
    sample_cnv_calls <- sample_cnv_calls %>%
      mutate(outlier_sample = ifelse(sample_id_trimmed %in% outliers, "YES", "NO"))
  } else {
    sample_cnv_calls <- sample_cnv_calls %>%
      mutate(outlier_sample = "NO")
  }
  
  # Filter based on exon type and select relevant columns
  sample_cnv_calls <- sample_cnv_calls %>%
    mutate(
      exon_type = case_when(
        Number.of.Exons.by.VarSeq.call == 1 ~ 'single',
        TRUE ~ 'multi'),
      # Remove 'chr' from VS_Call_chr
      VS_Call_chr = str_replace(VS_Call_chr, "^chr", ""),
      # Create a new 'coordinates' column by combining VS_Call_chr, VS_Call_start, and VS_Call_end
      coordinates = paste0(VS_Call_chr, ":", VS_Call_start, "-", VS_Call_end),
      log_p_value = round(-log(p.value.of.CNV, 10), 2),  # Log-transform and round
      Estimated.CN.of.CNV = factor(Estimated.CN.of.CNV, levels = c(0, 1, 2, 3, 4)),
      run = dnascreen_run
    )
  
  # Assign gene names
  sample_cnv_calls <- assign_gene_names(sample_cnv_calls, gene_data)
  
  # Check for overlaps with PLP CNVs
  sample_cnv_calls <- sample_cnv_calls %>%
    rowwise() %>%
    mutate(
      PLP = ifelse(any(sapply(1:nrow(plp_cnvs), function(i) {
        loc_split <- str_split(plp_cnvs$GRCh38Location[i], " - ")[[1]]
        start2 <- as.numeric(loc_split[1])
        end2 <- as.numeric(loc_split[2])
        check_overlap(
          chr1 = as.character(VS_Call_chr),
          start1 = VS_Call_start,
          end1 = VS_Call_end,
          chr2 = as.character(plp_cnvs$GRCh38Chromosome[i]),
          start2 = start2,
          end2 = end2
        )
      })), "YES", "NO")
    )
  
  # Adjust text for tooltip display to include Gene.Names and PLP, rounding numeric values
  sample_cnv_calls <- sample_cnv_calls %>%
    mutate(
      text = paste("Sample:", sample, "<br>",
                   "Coordinates:", coordinates, "<br>",
                   "Exon Type:", exon_type, "<br>",
                   "P-Value:", round(p.value.of.CNV, 4), "<br>",
                   "Estimated CN:", Estimated.CN.of.CNV, "<br>",
                   "Avg Z-Score of CNV:", round(Avg.Z.Score.of.CNV, 2), "<br>",
                   "Avg Read Ratio:", round(Avg.Read.Ratio.of.CNV, 2), "<br>",
                   "Run:", dnascreen_run, "<br>",
                   "Percent Difference:", round(Percent.Difference, 2), "<br>",
                   "Outlier Sample:", outlier_sample, "<br>",
                   "PLP:", PLP, "<br>",
                   "Gene Names:", ifelse(is.na(Gene.Names), "N/A", Gene.Names))
    )
  
  return(sample_cnv_calls)
}



#####################################################################################

# # Calculate the minimum -log p-values for CN=1 and CN=3 from simulated CNVs
# min_p_values_single <- sim_recalled_singleexon_table_combined %>%
#   filter(cn %in% c(1, 3)) %>%
#   summarise(
#     min_log_p_1 = min(log_p_value[cn == 1], na.rm = TRUE),
#     min_log_p_3 = min(log_p_value[cn == 3], na.rm = TRUE)
#   )
# 
# min_p_values_multi <- sim_recalled_clinvar_table_combined %>%
#   filter(cn %in% c(1, 3)) %>%
#   summarise(
#     min_log_p_1 = min(log_p_value[cn == 1], na.rm = TRUE),
#     min_log_p_3 = min(log_p_value[cn == 3], na.rm = TRUE)
#   )
# 
# # Extract the minimum values for easier reference
# min_log_p_del_single <- min_p_values_single$min_log_p_1
# min_log_p_dup_single <- min_p_values_single$min_log_p_3
# 
# min_log_p_del_multi <- min_p_values_multi$min_log_p_1
# min_log_p_dup_multi <- min_p_values_multi$min_log_p_3
# 
# # Create a color column in sample_cnv_calls_combined
# sample_cnv_calls_combined$color <- ifelse(sample_cnv_calls_combined$outlier_sample == "YES", "Outlier", "Normal")
# 
# # Create the single exon plot with jitter, horizontal lines, and annotations
# p_single <- ggplot() +
#   geom_jitter(data = sim_single, aes(x = cn, y = log_p_value, color = "Simulated CNVs"), alpha = 0.7, width = 0.2) +
#   geom_jitter(data = sample_cnv_calls_combined, aes(x = cn, y = log_adjusted_p_value, color = color), size = 2, alpha = 0.6, width = 0.2) +
#   geom_hline(yintercept = min_log_p_del_single, linetype = "dashed", color = "black", size = 1) +
#   geom_hline(yintercept = min_log_p_dup_single, linetype = "dashed", color = "black", size = 1) +
#   geom_text(aes(x = 1, y = min_log_p_del_single + 0.5, label = paste("CN=1:", round(min_log_p_del_single, 2))), 
#             vjust = -1, color = "black") +  # Adjust vertical positioning as needed
#   geom_text(aes(x = 3.5, y = min_log_p_dup_single + 0.5, label = paste("CN=3:", round(min_log_p_dup_single, 2))), 
#             vjust = 1, color = "black") +  # Adjust vertical positioning as needed
#   labs(
#     title = "Single Exon CNVs: -Log P-Values",
#     x = "Copy Number (CN)",
#     y = "-Log P-Values",
#     color = "Legend"
#   ) +
#   theme_minimal() +
#   scale_color_manual(values = c("Simulated CNVs" = "blue", "Normal" = "red", "Outlier" = "purple")) +
#   scale_x_discrete(breaks = c(0, 1, 2, 3, 4))
# 
# # Create the multi exon plot with jitter, horizontal lines, and annotations
# p_multi <- ggplot() +
#   geom_jitter(data = sim_multi, aes(x = cn, y = log_p_value, color = "Simulated CNVs"), alpha = 0.7, width = 0.2) +
#   geom_jitter(data = sample_cnv_calls_combined, aes(x = cn, y = log_adjusted_p_value, color = color), size = 2, alpha = 0.6, width = 0.2) +
#   geom_hline(yintercept = min_log_p_del_multi, linetype = "dashed", color = "black", size = 1) +
#   geom_hline(yintercept = min_log_p_dup_multi, linetype = "dashed", color = "black", size = 1) +
#   geom_text(aes(x = 1, y = min_log_p_del_multi + 0.5, label = paste("CN=1:", round(min_log_p_del_multi, 2))), 
#             vjust = -1, color = "black") +  # Adjust vertical positioning as needed
#   geom_text(aes(x = 3.5, y = min_log_p_dup_multi + 0.5, label = paste("CN=3:", round(min_log_p_dup_multi, 2))), 
#             vjust = 1, color = "black") +  # Adjust vertical positioning as needed
#   labs(
#     title = "Multi Exon CNVs: -Log P-Values",
#     x = "Copy Number (CN)",
#     y = "-Log P-Values",
#     color = "Legend"
#   ) +
#   theme_minimal() +
#   scale_color_manual(values = c("Simulated CNVs" = "blue", "Normal" = "red", "Outlier" = "purple")) +
#   scale_x_discrete(breaks = c(0, 1, 2, 3, 4))
# 
# # Display the plots
# p_single
# p_multi
# 
# 
# # Convert ggplot to ggplotly
# p_single_interactive <- ggplotly(p_single, tooltip = "text")
# p_multi_interactive <- ggplotly(p_multi, tooltip = "text")

# Save the interactive plot as an HTML file
# htmlwidgets::saveWidget(p_single_interactive, "summary_analyses_cnv_calls/p_single_interactive.html")
# htmlwidgets::saveWidget(p_multi_interactive, "summary_analyses_cnv_calls/p_multi_interactive.html")


# # Summarise p-values for simulation data
# p_values_df <- sim_table_combined %>%
#   group_by(cn, exon_type) %>%
#   summarise(
#     mean_p_value = mean(log_p_value, na.rm = TRUE),
#     ci_lower = mean(log_p_value, na.rm = TRUE) - qnorm(0.975) * sd(log_p_value, na.rm = FALSE) / sqrt(sum(!is.na(log_p_value))),
#     ci_upper = mean(log_p_value, na.rm = TRUE) + qnorm(0.975) * sd(log_p_value, na.rm = FALSE) / sqrt(sum(!is.na(log_p_value)))
#   )
#
# # Ensure cn includes all CN values as factors
# p_values_df <- p_values_df %>%
#   mutate(cn = factor(cn, levels = c(0, 1, 2, 3, 4)))


#
# # Define dodge width for alignment
# dodge_width <- 0.5
#
# # Plot the mean Z-Scores with 95% confidence intervals for each CN
# # along with CNVs called on real data of CNV calls made on real data
# p <- ggplot(p_values_df, aes(x = cn, y = mean_p_value, color = exon_type)) +
#   geom_point(aes(shape = "Simulated CNVs"),
#              position = position_dodge(width = dodge_width),
#              size = 3) +
#   geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
#                 width = 0.2, position = position_dodge(width = dodge_width)) +
#   geom_point(data = sample_cnv_calls_combined,
#              aes(x = cn, y = log_adjusted_p_value, color = exon_type,
#                  shape = "CNVs called on real data", text = text),
#              alpha = 0.6, size = 2,
#              position = position_dodge(width = dodge_width)) +
#   labs(
#     title = paste0("-Log adjusted p-values for Simulated CNV calls and CNVs called on real data"),
#     x = "Copy Number (CN)",
#     y = "-Log P-Values",
#     color = "Exon Type",
#     shape = "Legend"
#   ) +
#   theme_minimal() +
#   scale_color_manual(
#     values = c("single" = "blue", "multi" = "red"),
#     labels = c("single" = "Single Exon", "multi" = "Multi Exon")
#   ) +
#   scale_shape_manual(
#     values = c("Simulated CNVs" = 16, "CNVs called on real data" = 17)  # Use different shapes for differentiation
#   ) +
#   scale_x_discrete(breaks = c(0, 1, 3, 4)) +
#   theme(legend.position = "none")
#
# p
#
# # Convert ggplot to ggplotly
# p_interactive <- ggplotly(p, tooltip = "text")
# #
# # # Save the interactive plot as an HTML file
# # htmlwidgets::saveWidget(p_interactive,
# #                         "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/CI_plots/interactive_plot.html")
