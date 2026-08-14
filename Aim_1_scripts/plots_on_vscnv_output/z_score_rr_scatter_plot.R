## This plot is for Section 3 of my thesis where it outputs the 
## naive CNVs that are within the simulated range
## AND show which CNVs are selected for validation!

library(ggplot2)
library(plotly)
library(dplyr)
library(ggpubr)

real_data <- sample_cnv_calls_combined
# sent_for_acgh_cnvs <- read.csv("/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/Validations/files_sent_to_dan/cnvs_passed_three_cutoffs.csv")
# sent_for_acgh_cnvs$cn <- factor(sent_for_acgh_cnvs$cn, levels = c(0, 1, 2, 3, 4))
# 
# # Add validation column based on matching conditions
# real_data_check_validation <- real_data %>%
#   left_join(sent_for_acgh_cnvs, by = c("sample" = "sample_full_id", "cn", "coordinates")) %>%
#   mutate(sent_for_acgh = ifelse(is.na(exon), FALSE, TRUE)) %>%
#   select(-starts_with("extra_columns_from_acgh"))  # Remove unnecessary columns from sent_for_acgh_cnvs if needed
# 
# real_data$sent_for_acgh <- real_data_check_validation$sent_for_acgh
# # Since DNS001303 is requested for aCGH, set it TRUE for sent_for_acgh
# real_data[which(real_data$sample_id_trimmed == "DNS001303"), "sent_for_acgh"] <- TRUE

### Create a column in real_data to indicate whether CNV passed filter ###
# real_data <- real_data %>%
#   mutate(passed_score_filters = ifelse((cn == 3 & (log_adjusted_p_value > min_clinvar_dup_log_p_value) & 
#                                     (Avg.Read.Ratio.of.CNV > min_clinvar_dup_rr_diff) & 
#                                     (Avg.Z.Score.of.CNV > min_clinvar_dup_z_score)) |
#                                    (cn == 1 & (log_adjusted_p_value > min_clinvar_del_log_p_value) & 
#                                       (Avg.Read.Ratio.of.CNV < max_clinvar_del_rr_diff) & 
#                                       (Avg.Z.Score.of.CNV < max_clinvar_del_z_score)), TRUE, FALSE)) %>%
#   mutate(fail_percent_diff = ifelse(Percent.Difference > 19, TRUE, FALSE))

real_data <- real_data %>%
  mutate(passed_score_filters = ifelse((cn == 3 & 
                                          (Avg.Read.Ratio.of.CNV > min_clinvar_dup_rr_diff) & 
                                          (Avg.Z.Score.of.CNV > min_clinvar_dup_z_score)) |
                                         (cn == 1 & 
                                            (Avg.Read.Ratio.of.CNV < max_clinvar_del_rr_diff) & 
                                            (Avg.Z.Score.of.CNV < max_clinvar_del_z_score)), TRUE, FALSE)) %>%
  mutate(fail_percent_diff = ifelse(Percent.Difference > 19, TRUE, FALSE))


# filter out CNVs from outlier samples

real_data_without_outliers <- real_data %>%
  filter(outlier_sample == "NO") %>%
  filter(fail_percent_diff == FALSE)



real_data_without_outliers$color_group <- real_data_without_outliers$data_type



real_data_without_outliers <- real_data_without_outliers %>%
  group_by(data_type) %>%
  mutate(color_group = case_when(
    # sent_for_acgh == TRUE ~ "Selected for aCGH Validation",
    # fail_percent_diff == TRUE ~ "From sample with percent difference >20%",
    passed_score_filters == TRUE ~ "Naïve CNV calls that passed read ratio and z-score filters",
    color_group == "Real" ~ "Naïve CNV calls"
    # color_group != "False Positive" & !is.na(cn) ~ "Recalled simulated CNVs",
    # TRUE ~ "gray"
  ))

# #### Z SCORE vs p-value #### 
# # validation cnvs plot
# real_data_plot_z_score_vs_p_value <- ggplot(real_data_without_outliers,
#                                           aes(x = Avg.Z.Score.of.CNV,
#                                               y = log_adjusted_p_value,
#                                               color = color_group)) +
#   geom_point(size = 3, alpha = 0.6) +
#   
#   # # Vertical lines for Z-score minimums
#   # geom_vline(xintercept = c(max_clinvar_del_z_score, min_clinvar_dup_z_score),
#   #            linetype = "dashed", color = "black") +
#   # 
#   # # Partial horizontal lines for deletions and duplications (colored red)
#   # annotate("segment", x = -15, xend = 0, y = min_clinvar_del_log_p_value,
#   #          yend = min_clinvar_del_log_p_value, linetype = "dashed", color = "red") +
#   # annotate("segment", x = 0, xend = 15, y = min_clinvar_dup_log_p_value,
#   #          yend = min_clinvar_dup_log_p_value, linetype = "dashed", color = "red") +
#   # 
#   # # Compact annotations for minimum values
#   # annotate("text", x = max_clinvar_del_z_score - 7,
#   #          y = max(real_data_without_outliers$log_adjusted_p_value),
#   #          label = paste0("deletion z-score maximum = ", round(max_clinvar_del_z_score, 4)),
#   #          color = "black", vjust = -1) +
#   # annotate("text", x = min_clinvar_dup_z_score + 7,
#   #          y = max(real_data_without_outliers$log_adjusted_p_value),
#   #          label = paste0("duplication z-score minimum = ", round(min_clinvar_dup_z_score, 4)),
#   #          color = "black", vjust = -1) +
#   # annotate("text", x = -12,
#   #          y = min_clinvar_del_log_p_value,
#   #          label = paste0("deletion log p minimum = ", round(min_clinvar_del_log_p_value, 4)),
#   #          color = "red", vjust = 1.5) +
#   # annotate("text", x = 12,
#   #          y = min_clinvar_dup_log_p_value,
#   #          label = paste0("duplication log p minimum = ", round(min_clinvar_dup_log_p_value, 4)),
#   #          color = "red", vjust = -1.5) +
#   
#   # Labels and styling
#   labs(
#     # title = "Z-Score vs -log(p-value) of Naïve CNV Calls",
#     x = "Z-Score",
#     y = "-log(p-value)",
#     color = "CNV Type: "
#   ) +
#   scale_color_manual(values = c(
#     'Naïve CNV calls that passed p-value, read ratio and z-score filters' = 'pink',
#     'Naïve CNV calls' = 'pink',
#     'Selected for aCGH Validation' = 'pink',
#     'Recalled simulated CNVs' = '#A6C8FF'
#   )) +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "none",
#     legend.direction = "horizontal",
#     legend.title = element_text(size = 13, face = "bold"),
#     legend.text = element_text(size = 12),
#     legend.key.size = unit(0.8, "lines"),
#     legend.spacing.x = unit(0.4, 'cm'),
#     
#     plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#     axis.title.x = element_text(size = 14, face = "bold"),
#     axis.title.y = element_text(size = 14, face = "bold"),
#     axis.text = element_text(size = 12)
#   )
# 
# real_data_plot_z_score_vs_p_value
# 
# # Convert to ggplotly
# # real_data_plot_rr_vs_p_value_plotly <- ggplotly(real_data_plot_rr_vs_p_value, tooltip = "text")
# 
# 
# 
# #### READ RATIO vs p-value #### 
# # validation cnvs plot
# real_data_plot_rr_vs_p_value <- ggplot(real_data_without_outliers,
#                                           aes(x = Avg.Read.Ratio.of.CNV,
#                                               y = log_adjusted_p_value,
#                                               color = color_group)) +
#   geom_point(size = 3, alpha = 0.6) +
#   
#   # # Vertical lines for max and min Read Ratios
#   # geom_vline(xintercept = max_clinvar_del_rr_diff, linetype = "dashed", color = "white") +
#   # geom_vline(xintercept = min_clinvar_dup_rr_diff, linetype = "dashed", color = "white") +
#   # 
#   # # Partial horizontal lines for deletions and duplications (colored red)
#   # annotate("segment", x = 1, xend = 0, y = min_clinvar_del_log_p_value,
#   #          yend = min_clinvar_del_log_p_value, linetype = "dashed", color = "white") +
#   # annotate("segment", x = 1, xend = 2, y = min_clinvar_dup_log_p_value,
#   #          yend = min_clinvar_dup_log_p_value, linetype = "dashed", color = "white") +
# 
#   # Annotations for minimums, positioned further from the center
#   # annotate("text", x = max_clinvar_del_rr_diff - 0.5,
#   #          y = max(real_data_without_outliers$log_p_value),
#   #          label = paste0("deletion read ratio maximum = ", round(max_clinvar_del_rr_diff, 4)),
#   #          color = "black", vjust = -1) +
#   # annotate("text", x = min_clinvar_dup_rr_diff + 0.5,
#   #          y = max(real_data_without_outliers$log_p_value),
#   #          label = paste0("duplication read ratio minimum = ", round(min_clinvar_dup_rr_diff, 4)),
#   #          color = "black", vjust = -1) +
#   # annotate("text", x = 0.5,
#   #          y = min_clinvar_del_log_p_value,
#   #          label = paste0("deletion log p minimum = ", round(min_clinvar_del_log_p_value, 4)),
#   #          color = "red", vjust = 2) +
#   # annotate("text", x = 1.5,
#   #          y = min_clinvar_dup_log_p_value,
#   #          label = paste0("duplication log p minimum = ", round(min_clinvar_dup_log_p_value, 4)),
#   #          color = "red", vjust = 2) +
#   
#   # Labels and styling
#   labs(
#     # title = "Read Ratio vs -log(p-value) of Naïve CNV Calls",
#     x = "Read Ratio",
#     y = "-log(p-value)",
#     color = "CNV Type: "
#   ) +
#   scale_color_manual(values = c(
#     'Naïve CNV calls that passed p-value, read ratio and z-score filters' = 'pink',
#     'Naïve CNV calls' = 'pink',
#     'Selected for aCGH Validation' = 'pink',
#     'Recalled simulated CNVs' = '#A6C8FF'
#   )) +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "none",
#     legend.direction = "horizontal",
#     legend.title = element_text(size = 13, face = "bold"),
#     legend.text = element_text(size = 12),
#     legend.key.size = unit(0.8, "lines"),
#     legend.spacing.x = unit(0.4, 'cm'),
#     
#     plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#     axis.title.x = element_text(size = 14, face = "bold"),
#     axis.title.y = element_text(size = 14, face = "bold"),
#     axis.text = element_text(size = 12)
#   )
# 
# real_data_plot_rr_vs_p_value

# Convert to ggplotly
# real_data_plot_rr_vs_p_value_plotly <- ggplotly(real_data_plot_rr_vs_p_value, tooltip = "text")


#### Calculate deviation counts/percentages ####

deletion_deviation_stats <- real_data_without_outliers %>%
  filter(Avg.Read.Ratio.of.CNV < 1.0) %>%
  summarise(
    total = n(),
    outside_tolerance = sum(Avg.Read.Ratio.of.CNV < 0.4 | Avg.Read.Ratio.of.CNV > 0.6, na.rm = TRUE),
    pct_outside = round(100 * outside_tolerance / total, 1)
  )

duplication_deviation_stats <- real_data_without_outliers %>%
  filter(Avg.Read.Ratio.of.CNV > 1.0) %>%
  summarise(
    total = n(),
    outside_tolerance = sum(Avg.Read.Ratio.of.CNV < 1.4 | Avg.Read.Ratio.of.CNV > 1.6, na.rm = TRUE),
    pct_outside = round(100 * outside_tolerance / total, 1)
  )

print(paste0("Deletions outside ±0.1 of 0.5: ", deletion_deviation_stats$outside_tolerance,
             " of ", deletion_deviation_stats$total,
             " (", deletion_deviation_stats$pct_outside, "%)"))

print(paste0("Duplications outside ±0.1 of 1.5: ", duplication_deviation_stats$outside_tolerance,
             " of ", duplication_deviation_stats$total,
             " (", duplication_deviation_stats$pct_outside, "%)"))





#### READ RATIO vs Z-SCORE ####

#### Plot ####

real_data_plot_rr_vs_z_score <- ggplot(real_data_without_outliers,
                                       aes(x = Avg.Read.Ratio.of.CNV,
                                           y = Avg.Z.Score.of.CNV,
                                           color = color_group, text = text)) +
  geom_point(size = 3, alpha = 0.6) +
  
  # Thicker reference lines at Z-score = 0 and Read Ratio = 1.0
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.25) +
  geom_vline(xintercept = 1.0, color = "grey40", linewidth = 0.25) +
  
  # Theoretical deletion/duplication read ratio lines (blue dotted)
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "blue", linewidth = 0.6) +
  geom_vline(xintercept = 1.5, linetype = "dotted", color = "blue", linewidth = 0.6) +
  
  # ±0.1 tolerance bands (light blue dotted)
  geom_vline(xintercept = 0.4, linetype = "dotted", color = "blue", linewidth = 0.5) +
  geom_vline(xintercept = 0.6, linetype = "dotted", color = "blue", linewidth = 0.5) +
  geom_vline(xintercept = 1.4, linetype = "dotted", color = "blue", linewidth = 0.5) +
  geom_vline(xintercept = 1.6, linetype = "dotted", color = "blue", linewidth = 0.5) +
  
  # Leader lines connecting text to the theoretical dotted lines
  annotate("segment",
           x = 0.5 - 0.03, xend = 0.5,
           y = max(real_data_without_outliers$Avg.Z.Score.of.CNV, na.rm = TRUE),
           yend = max(real_data_without_outliers$Avg.Z.Score.of.CNV, na.rm = TRUE)-1,
           color = "blue", linewidth = 0.4) +
  annotate("segment",
           x = 1.5 + 0.03, xend = 1.5,
           y = max(real_data_without_outliers$Avg.Z.Score.of.CNV, na.rm = TRUE),
           yend = max(real_data_without_outliers$Avg.Z.Score.of.CNV, na.rm = TRUE)-1,
           color = "blue", linewidth = 0.4) +
  
  # Annotations for theoretical lines (horizontal text, positioned to the side)
  annotate("text", x = 0.5 - 0.03,
           y = max(real_data_without_outliers$Avg.Z.Score.of.CNV, na.rm = TRUE),
           label = "Theoretical deletion ratio = 0.5", color = "blue",
           hjust = 1, vjust = 0.5, size = 4) +
  annotate("text", x = 1.5 + 0.03,
           y = max(real_data_without_outliers$Avg.Z.Score.of.CNV, na.rm = TRUE),
           label = "Theoretical duplication ratio = 1.5", color = "blue",
           hjust = 0, vjust = 0.5, size = 4) +
  
  # Custom x-axis breaks
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8)) +
  
  # Labels and styling
  labs(
    x = "Read Ratio",
    y = "Z-Score",
    color = "CNV Type: "
  ) +
  scale_color_manual(values = c(
    'Naïve CNV calls that passed read ratio and z-score filters' = 'pink',
    'Naïve CNV calls' = 'pink',
    'Selected for aCGH Validation' = 'pink'
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
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

real_data_plot_rr_vs_z_score

# Convert to ggplotly
# real_data_plot_rr_vs_z_score_plotly <- ggplotly(real_data_plot_rr_vs_z_score, tooltip = "text")
# 
# real_data_plot_rr_vs_z_score_plotly




