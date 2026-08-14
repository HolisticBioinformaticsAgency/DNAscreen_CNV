## This plot is for Section 3 of my thesis where it outputs the 
## VS-CNV calls that are within the simulated range.

###### RUN filtering_naive_calls_for_validation.R FIRST! ####

library(ggplot2)
library(plotly)
library(dplyr)
library(ggpubr)

sim_clinvar_del <- sim_recalled_clinvar_table_combined %>%
  filter(cn == 0 | cn == 1)

sim_clinvar_dup <- sim_recalled_clinvar_table_combined %>%
  filter(cn == 3 | cn == 4)

max_clinvar_del_z_score <- max(sim_clinvar_del$simulated_z_score)
min_clinvar_dup_z_score <- min(sim_clinvar_dup$simulated_z_score)

max_clinvar_del_rr_diff <- max(sim_clinvar_del$Avg.Read.Ratio.of.CNV)
min_clinvar_dup_rr_diff <- min(sim_clinvar_dup$Avg.Read.Ratio.of.CNV)

min_clinvar_del_log_p_value <- min(sim_clinvar_del$log_p_value)
min_clinvar_dup_log_p_value <- min(sim_clinvar_dup$log_p_value)

### Create a column in real_data to indicate whether CNV passed filter ###
real_data <- real_data %>%
  mutate(fail_percent_diff = ifelse(Percent.Difference >= 20, TRUE, FALSE))


# for the clinvar, we want to remove the CNVs that are not reliably detected by VS-CNV
# collapsed_sim_missed_clinvar_table_combined derived from clinvar_cnvs_exon_percentage_overlap_analysis.R
sim_recalled_clinvar_table_combined <- anti_join(sim_recalled_clinvar_table_combined, 
                                                 collapsed_sim_missed_clinvar_table_combined, 
                                                 by = c("chr", "start", "end", "gene"))

sim_recalled_clinvar_table_combined <- sim_recalled_clinvar_table_combined %>%
  mutate(Avg.Z.Score.of.CNV = simulated_z_score)

# filter out CNVs from samples > 20 avg percent diff as they
# have many false positives
real_data_without_outliers <- real_data %>%
  filter(outlier_sample == "NO")

combined_clinvar_data <- bind_rows(sim_recalled_clinvar_table_combined,
                                   real_data_without_outliers)


combined_clinvar_data$color_group <- combined_clinvar_data$data_type


### Z SCORE vs P-value ###


combined_clinvar_data <- combined_clinvar_data %>%
  group_by(data_type) %>%
  mutate(color_group = case_when(
    data_type == "Simulated" ~ "Recalled simulated CNVs",
    # (sent_for_acgh == TRUE & passed_score_filters == TRUE) ~ "Passed and sent for aCGH validation",
    # (sent_for_acgh == TRUE & passed_score_filters == FALSE) ~ "Failed and sent for aCGH validation",
    passed_z_score == TRUE & passed_pval == TRUE ~ "VS-CNV calls that are within Z-Score and p-value range of simulated CNVs",
    data_type == "Real" ~ "VS-CNV calls",
    TRUE ~ "gray"
  ))

combined_clinvar_data <- combined_clinvar_data %>%
  mutate(color_group = factor(color_group))

# Clinvar CNV plot with colored dots by copy number and false positives in red
clinvar_exon_plot_z_score_vs_p_value <- ggplot(combined_clinvar_data,
                                               aes(x = Avg.Z.Score.of.CNV,
                                                   y = log_adjusted_p_value,
                                                   color = color_group,
                                                   text = text)) +
  geom_point(size = 3, alpha = 0.6) +

  # Vertical lines for Z-score minimums
  geom_vline(xintercept = c(max_clinvar_del_z_score, min_clinvar_dup_z_score),
             linetype = "dashed", color = "black") +

  # Partial horizontal lines for deletions and duplications (colored red)
  annotate("segment", x = -23, xend = 0, y = min_clinvar_del_log_p_value,
           yend = min_clinvar_del_log_p_value, linetype = "dashed", color = "red") +
  annotate("segment", x = 0, xend = 15, y = min_clinvar_dup_log_p_value,
           yend = min_clinvar_dup_log_p_value, linetype = "dashed", color = "red") +

  # Compact annotations for minimum values
  annotate("text", x = max_clinvar_del_z_score - 7,
           y = max(combined_clinvar_data$log_adjusted_p_value),
           label = paste0("deletion z-score maximum = ", round(max_clinvar_del_z_score, 2)),
           color = "black", vjust = -1) +
  annotate("text", x = min_clinvar_dup_z_score + 7,
           y = max(combined_clinvar_data$log_adjusted_p_value),
           label = paste0("duplication z-score minimum = ", round(min_clinvar_dup_z_score, 2)),
           color = "black", vjust = -1) +
  annotate("text", x = -12,
           y = min_clinvar_del_log_p_value,
           label = paste0("deletion log p minimum = ", round(min_clinvar_del_log_p_value, 2)),
           color = "red", vjust = 1.5) +
  annotate("text", x = 12,
           y = min_clinvar_dup_log_p_value,
           label = paste0("duplication log p minimum = ", round(min_clinvar_dup_log_p_value, 2)),
           color = "red", vjust = -1.5) +

  # Labels and styling
  labs(
    title = "Comparison of Z-Score vs -log(p-value) for Simulated CNVs and VS-CNV calls",
    x = "Z-Score",
    y = "-log(p-value)",
    color = "CNV Type: "
  ) +
  scale_color_manual(values = c(
    # 'Failed and sent for aCGH validation' = 'pink',
    'VS-CNV calls' = 'pink',
    'VS-CNV calls that are within Z-Score and p-value range of simulated CNVs' = 'green',
    # 'Failed confidence metric filters' = 'pink',
    # 'Passed and sent for aCGH validation' = 'green',
    'Recalled simulated CNVs' = '#A6C8FF'
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.8, "lines"),
    legend.spacing.x = unit(0.4, 'cm'),
    
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )


# Display the plot
clinvar_exon_plot_z_score_vs_p_value
# 
# 
# # Convert to ggplotly
# clinvar_exon_plot_z_score_vs_p_value_plotly <- ggplotly(clinvar_exon_plot_z_score_vs_p_value, tooltip = "text")


#### READ RATIO vs P-value ####


combined_clinvar_data <- combined_clinvar_data %>%
  group_by(data_type) %>%
  mutate(color_group = case_when(
    data_type == "Simulated" ~ "Recalled simulated CNVs",
    # (sent_for_acgh == TRUE & passed_score_filters == TRUE) ~ "Passed and sent for aCGH validation",
    # (sent_for_acgh == TRUE & passed_score_filters == FALSE) ~ "Failed and sent for aCGH validation",
    passed_rr == TRUE & passed_pval == TRUE ~ "VS-CNV calls that are within read ratio and p-value range of simulated CNVs",
    data_type == "Real" ~ "VS-CNV calls",
    TRUE ~ "gray"
  ))

combined_clinvar_data <- combined_clinvar_data %>%
  mutate(color_group = factor(color_group))
  

# clinvar exon plot
clinvar_exon_plot_rr_vs_p_value <- ggplot() +
  geom_point(data = combined_clinvar_data,
             aes(x = Avg.Read.Ratio.of.CNV,
                 y = log_adjusted_p_value,
                 color = color_group,
                 text = text),
             size = 3, alpha = 0.6) +
  
  # Vertical lines for max and min Read Ratios
  geom_vline(xintercept = max_clinvar_del_rr_diff, linetype = "dashed", color = "black") +
  geom_vline(xintercept = min_clinvar_dup_rr_diff, linetype = "dashed", color = "black") +
  
  # Partial horizontal lines for deletions and duplications (colored red)
  annotate("segment", x = 1, xend = 0, y = min_clinvar_del_log_p_value,
           yend = min_clinvar_del_log_p_value, linetype = "dashed", color = "black") +
  annotate("segment", x = 1, xend = 2, y = min_clinvar_dup_log_p_value,
           yend = min_clinvar_dup_log_p_value, linetype = "dashed", color = "black") +

  # Annotations for minimums, positioned further from the center
  annotate("text", x = max_clinvar_del_rr_diff - 0.3,
           y = max(combined_clinvar_data$log_p_value),
           label = paste0("deletion read ratio maximum = ", round(max_clinvar_del_rr_diff, 2)),
           color = "black", vjust = -1) +
  annotate("text", x = min_clinvar_dup_rr_diff + 0.3,
           y = max(combined_clinvar_data$log_p_value),
           label = paste0("duplication read ratio minimum = ", round(min_clinvar_dup_rr_diff, 2)),
           color = "black", vjust = -1) +
  annotate("text", x = 0.5,
           y = min_clinvar_del_log_p_value,
           label = paste0("deletion log p minimum = ", round(min_clinvar_del_log_p_value, 2)),
           color = "red", vjust = 2) +
  annotate("text", x = 1.5,
           y = min_clinvar_dup_log_p_value,
           label = paste0("duplication log p minimum = ", round(min_clinvar_dup_log_p_value, 2)),
           color = "red", vjust = 2) +

  # Labels and styling
  labs(
    title = "Comparison of Read Ratio vs -log(p-value) for Simulated CNVs and VS-CNV calls",
    x = "Read Ratio",
    y = "-log(p-value)",
    color = "CNV Type: ",
    size = 22
  ) +
  scale_color_manual(values = c(
    # 'Failed and sent for aCGH validation' = 'pink',
    'VS-CNV calls' = 'pink',
    'VS-CNV calls that are within read ratio and p-value range of simulated CNVs' = 'green',
    # 'Failed confidence metric filters' = 'pink',
    # 'Passed and sent for aCGH validation' = 'green',
    'Recalled simulated CNVs' = '#A6C8FF'
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.8, "lines"),
    legend.spacing.x = unit(0.4, 'cm'),
    
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )


clinvar_exon_plot_rr_vs_p_value

# Convert to ggplotly
# clinvar_exon_plot_rr_vs_p_value_plotly <- ggplotly(clinvar_exon_plot_rr_vs_p_value, tooltip = "text")



