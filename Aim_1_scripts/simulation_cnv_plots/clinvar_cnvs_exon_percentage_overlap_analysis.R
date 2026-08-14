library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)

# Investigating CNVs encompassing one exon
sim_recalled_clinvar_encompassing_one_exon <- sim_recalled_clinvar_table_combined %>%
  filter(count_exons_within_simulated_window == 1)

# Investigating CNVs encompassing more than one exon
sim_recalled_clinvar_encompassing_more_than_one_exon <- sim_recalled_clinvar_table_combined %>%
  filter(count_exons_within_simulated_window > 1)

# Investigating only recalled CNVs that have no exon (capture target region) fully within the coordinates of the simulated CNV
# to compare with the missed CNVs (ALL MISSED CNVs have no exon fully within coordinates)
# Filter recalled CNVs with no exon fully within coordinates
sim_recalled_clinvar_no_exon_within <- sim_recalled_clinvar_table_combined %>%
  filter(count_exons_within_simulated_window == 0)

### Finding difference in terms of Read ratio changes from missed vs recalled CNVs
# sim_clinvar_partially_overlapping_table_combined <- bind_rows(sim_missed_clinvar_table_combined,
#                                                               sim_recalled_clinvar_no_exon_within)


rr_z_diff_missed <- sim_missed_clinvar_table_combined %>%
  mutate(
    rr_difference = Read.Ratio.of.Simulated.CNV - orig_rr,
    z_difference  = simulated_z_score - orig_z_score,
    overlap_class = "Missed CNVs "
  )

# Recalled, partially overlapping CNVs (0 exons within)
rr_z_diff_partial <- sim_recalled_clinvar_table_combined %>%
  mutate(
    rr_difference = Read.Ratio.of.Simulated.CNV - orig_rr,
    z_difference  = simulated_z_score - orig_z_score,
    overlap_class = "Recalled CNVs"
  )

rr_z_diff_combined <- bind_rows(
  rr_z_diff_missed,
  rr_z_diff_partial
)

### histogram of Δ read ratio of missed vs recalled CNVs

# calculate bin width once (keeps scaling consistent)
binwidth_val <- diff(range(rr_z_diff_combined$rr_difference, na.rm = TRUE)) / 50
n_vals <- nrow(rr_z_diff_combined)

ggplot(rr_z_diff_combined, aes(x = rr_difference)) +
  
  # Histogram (counts)
  geom_histogram(
    aes(fill = overlap_class),
    bins = 50,
    alpha = 0.4,
    position = "identity"
  ) +
  
  # Density line (scaled to counts)
  geom_density(
    aes(
      y = after_stat(density * n_vals * binwidth_val),
      colour = overlap_class
    ),
    linewidth = 1.2
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  scale_y_continuous(
    name = "Number of CNVs",
    sec.axis = sec_axis(
      ~ . / (n_vals * binwidth_val),
      name = "Density"
    )
  ) +
  
  labs(
    x = "Δ Read Ratio (Post - Pre simulation)",
    fill = "CNV category",
    colour = "CNV category",
    title = "Δ Read Ratio (Post - Pre Simulation) of\nMissed vs Recalled CNVs"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "bottom"
  )

### plot of Δ read ratio vs exon overlap percentage
### Also shows minimum overlap
min_overlap_recalled <- rr_z_diff_combined %>%
  filter(overlap_class == "Recalled CNVs") %>%
  summarise(min_overlap = min(exon_overlap_percentage, na.rm = TRUE)) %>%
  pull(min_overlap)

ggplot(
  rr_z_diff_combined,
  aes(
    x = exon_overlap_percentage,
    y = rr_difference,
    colour = overlap_class
  )
) +
  geom_point(alpha = 0.6, size = 1.5) +
  
  # Horizontal line at y=0
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # Vertical line at minimum recalled exon overlap
  geom_vline(xintercept = min_overlap_recalled, linetype = "dotted", colour = "black", size = 0.8) +
  
  labs(
    x = "Percentage of Exonic Overlap",
    y = "Δ Read Ratio (Post - Pre Simulation)",
    colour = "CNV category",
    title = "Δ Read Ratio (Post - Pre Simulation) vs Percentage of Exonic Overlap\nMissed vs Recalled CNVs"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "bottom"
  ) + annotate(
    "text",
    x = min_overlap_recalled + 2, # shift slightly right
    y = max(rr_z_diff_combined$rr_difference, na.rm = TRUE),
    label = paste0("Minumum overlap for recalled CNVs = ", round(min_overlap_recalled, 1), "%"),
    colour = "black",
    hjust = 0
  )



### Below is to find minimum overlap where partially overlapped CNVs
### are consistently missed.

# Collapse rows and calculate frequency
collapse_cnvs_table <- function(input_df, missed=FALSE) {
  input_df %>%
    group_by(
      chr, start, end, gene, 
      overlap_with_start_name,
      overlap_with_end_name,
      cn,
      summed_length_exons_overlapping,
      within_exon,
      count_exons_overlapping_simulated_window,
      exon_overlap_length,
      exon_overlap_percentage
    ) %>%
    summarise(
      num_of_recalled_runs = ifelse(missed, 6 - n(), n()), # Count the frequency of each unique combination
      .groups = "drop" # Avoid creating grouped data in the output
    ) %>%
    mutate(
      exon_overlap_percentage = round(exon_overlap_percentage, 2),
      size_of_simulated_CNV = end - start,
      text = paste(
        "Coordinates:", paste0(chr, ": ", start, '-', end), "<br>",
        "Size of Simulated CNV:", size_of_simulated_CNV, "<br>",
        "Total overlap length with overlapping exons:", round(exon_overlap_length, 2), "<br>",
        "Percentage total length of overlapping exons that overlaps with simulated window:",
        exon_overlap_percentage, "<br>",
        "Number of Runs recalled in:", num_of_recalled_runs, "<br>",
        "Gene:", gene, "<br>",
        "Number of exons overlapping simulated window:", count_exons_overlapping_simulated_window, "<br>",
        "Left-hand Exon:", overlap_with_start_name, "<br>",
        "Right-hand Exon:", overlap_with_end_name, "<br>",
        "Summed Length of overlapping exons:", summed_length_exons_overlapping, "<br>",
        "Simulated CN:", cn, "<br>",
        "Within Exon:", within_exon, "<br>"
      )
    )
}

# CNVs that encompass one exon or more
# We are looking at recall table because i checked that there are no missed CNVs
# that encompass one exon or more!
collapsed_sim_recalled_clinvar_encompassing_one_exon <- 
  collapse_cnvs_table(sim_recalled_clinvar_encompassing_one_exon)
# write.csv(collapsed_sim_recalled_clinvar_encompassing_one_exon, 'summary_analyses_cnv_calls/simulation_summary_analyses/collapsed_sim_recalled_clinvar_encompassing_one_exon.csv')

collapsed_sim_recalled_clinvar_encompassing_more_than_one_exon <- 
  collapse_cnvs_table(sim_recalled_clinvar_encompassing_more_than_one_exon)
# write.csv(collapsed_sim_recalled_clinvar_encompassing_more_than_one_exon, 'summary_analyses_cnv_calls/simulation_summary_analyses/collapsed_sim_recalled_clinvar_encompassing_more_than_one_exon.csv')

# CNVs that DO NOT encompass an exon
collapsed_sim_recalled_clinvar_no_exon_within <- 
  collapse_cnvs_table(sim_recalled_clinvar_no_exon_within)

collapsed_sim_missed_clinvar_table_combined <- 
  collapse_cnvs_table(sim_missed_clinvar_table_combined, TRUE)


# Since CNVs that are missed are all CNVs that DO NOT encompass an exon
# The following analysis looks at only at CNVs of this category.

# missed_in_all <- collapsed_sim_missed_clinvar_table_combined %>%
#   filter(num_of_recalled_runs == 0) %>%
#   mutate(type = 'Missed in all six runs')
# collapsed_sim_recalled_clinvar_no_exon_within <- 
#   collapsed_sim_recalled_clinvar_no_exon_within %>%
#   mutate(type = 'Recalled in at least one run')

cnvs_not_encompass_exon <- bind_rows(
  collapsed_sim_recalled_clinvar_no_exon_within %>% mutate(source = "recalled"),
  collapsed_sim_missed_clinvar_table_combined %>% mutate(source = "missed")
) %>%
  arrange(source) %>%  # Ensures 'recalled' comes first
  distinct(chr, start, end, gene, cn, .keep_all = TRUE) %>%  # Keeps first occurrence
  select(-source)  # Removes helper column
# write.csv(cnvs_not_encompass_exon, "summary_analyses_cnv_calls/simulation_summary_analyses/cnvs_not_encompass_exon.csv")

cnvs_partially_overlap_one_exon <- cnvs_not_encompass_exon %>%
  filter(count_exons_overlapping_simulated_window == 1 & within_exon == FALSE)
# write.csv(cnvs_partially_overlap_one_exon, "summary_analyses_cnv_calls/simulation_summary_analyses/cnvs_partially_overlap_one_exon.csv")

cnvs_partially_overlap_two_exon <- cnvs_not_encompass_exon %>%
  filter(count_exons_overlapping_simulated_window > 1 & within_exon == FALSE)
# write.csv(cnvs_partially_overlap_two_exon, "summary_analyses_cnv_calls/simulation_summary_analyses/cnvs_partially_overlap_two_exon.csv")

cnvs_within_exon <- cnvs_not_encompass_exon %>%
  filter(within_exon == TRUE)
# write.csv(cnvs_within_exon, "summary_analyses_cnv_calls/simulation_summary_analyses/cnvs_within_exon.csv")


# Calculate the minimum exon_overlap_percentage for the specified subsets
min_recalled_partially_overlap_one_exon <- round(min(
  cnvs_partially_overlap_one_exon$exon_overlap_percentage
  [cnvs_partially_overlap_one_exon$num_of_recalled_runs > 0]), 2)

min_recalled_partially_overlap_more_than_one_exon <- round(min(
  cnvs_partially_overlap_two_exon$exon_overlap_percentage
  [cnvs_partially_overlap_two_exon$num_of_recalled_runs > 0]), 2)
                                  
min_recalled_within_exon <- round(min(
  cnvs_within_exon$exon_overlap_percentage
  [cnvs_within_exon$num_of_recalled_runs > 0]), 2)


# dot plot for each category above
#############################################
dot_plot_cnvs_partially_overlap_one_exon <- ggplot(cnvs_partially_overlap_one_exon, aes(
  x = exon_overlap_percentage, 
  y = exon_overlap_length, 
  color = as.factor(num_of_recalled_runs), 
  text = text
)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_y_log10() +  # Apply log scale to y-axis
  labs(
    title = "Percentage length of exon/s that are overlapping with the CNV and total number of exonic bases within the CNV
    Category: CNV overlapping ONE intron-exon boundary",
    x = "Percentage length of exon/s that are overlapping with the CNV",
    y = "Total number of exonic bases (log scale)",
    color = "Number of Recalled Runs"
  ) +
  theme_minimal() +
  scale_color_manual(
    values = c(
      "0" = "red",           # Red for 0
      "1" = "#dbe9f6",       # Light blue
      "2" = "#b4d3eb",       # Slightly darker blue
      "3" = "#8cbee1",       # Medium blue
      "4" = "#6699d6",       # Darker blue
      "5" = "#3f74cc",       # Even darker blue
      "6" = "#174fc1"        # Darkest blue
    )
  ) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +  # Circle for FALSE, Triangle for TRUE
  geom_vline(aes(
    xintercept = min_recalled_partially_overlap_one_exon,
    text = paste0("Minimum percentage overlap of recalled CNVs = ", min_recalled_partially_overlap_one_exon)
  ), color = "#1f77b4", linetype = "dashed", size = 0.8, alpha = 0.7) +
  annotate("text", 
           x = min_recalled_partially_overlap_one_exon + 1, 
           y = max(cnvs_partially_overlap_one_exon$exon_overlap_length, na.rm = TRUE) / 2,
           label = paste0("Minimum percentage overlap of recalled CNVs = ", min_recalled_partially_overlap_one_exon),
           color = "#1f77b4", hjust = 0, size = 4)

dot_plot_cnvs_partially_overlap_one_exon

#############################################

dot_plot_cnvs_partially_overlap_two_exon <- ggplot(cnvs_partially_overlap_two_exon, aes(x = exon_overlap_percentage, y = exon_overlap_length, 
                                      color = as.factor(num_of_recalled_runs), text = text)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_y_log10() +  # Apply log scale to y-axis
  labs(
    title = "Percentage length of exon/s that are overlapping with the CNV and total number of exonic bases within the CNV
    Category: CNV overlapping TWO intron-exon boundary",
    x = "Percentage length of exon/s that are overlapping with the CNV",
    y = "Total number of exonic bases (log scale)",
    color = "Number of Recalled Runs"
  ) +
  theme_minimal() +
  scale_color_manual(
    values = c(
      "0" = "red",           # Red for 0
      "1" = "#dbe9f6",       # Light blue
      "2" = "#b4d3eb",       # Slightly darker blue
      "3" = "#8cbee1",       # Medium blue
      "4" = "#6699d6",       # Darker blue
      "5" = "#3f74cc",       # Even darker blue
      "6" = "#174fc1"        # Darkest blue
    )
  ) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +  # Circle for FALSE, Triangle for TRUE
  geom_vline(aes(xintercept = min_recalled_partially_overlap_more_than_one_exon,
                 text = paste0("Minimum percentage overlap of recalled CNVs = ", min_recalled_partially_overlap_more_than_one_exon)),
             color = "#1f77b4", linetype = "dashed", size = 0.8, alpha = 0.7) +
  annotate("text", x = min_recalled_partially_overlap_more_than_one_exon + 1, y = max(cnvs_partially_overlap_two_exon$exon_overlap_length, na.rm = TRUE) / 2,
           label = paste0("Minimum percentage overlap of recalled CNVs = ", min_recalled_partially_overlap_more_than_one_exon),
           color = "#1f77b4", hjust = 0, size = 4)

dot_plot_cnvs_partially_overlap_two_exon


#############################################

dot_plot_cnvs_within_exon <- ggplot(cnvs_within_exon, aes(x = exon_overlap_percentage, y = exon_overlap_length, 
                                                          color = as.factor(num_of_recalled_runs), text = text)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_y_log10() +  # Apply log scale to y-axis
  labs(
    title = "Percentage length of exon/s that are overlapping with the CNV and total number of exonic bases within the CNV
    Category: Intra-exon CNV",
    x = "Percentage length of exon/s that are overlapping with the CNV",
    y = "Total number of exonic bases (log scale)",
    color = "Number of Recalled Runs"
  ) +
  theme_minimal() +
  scale_color_manual(
    values = c(
      "0" = "red",           # Red for 0
      "1" = "#dbe9f6",       # Light blue
      "2" = "#b4d3eb",       # Slightly darker blue
      "3" = "#8cbee1",       # Medium blue
      "4" = "#6699d6",       # Darker blue
      "5" = "#3f74cc",       # Even darker blue
      "6" = "#174fc1"        # Darkest blue
    )
  ) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +  # Circle for FALSE, Triangle for TRUE
  geom_vline(aes(xintercept = min_recalled_within_exon,
                 text = paste0("Minimum percentage overlap of recalled CNVs = ", min_recalled_within_exon)),
             color = "#1f77b4", linetype = "dashed", size = 0.8, alpha = 0.7) +
  annotate("text", x = min_recalled_within_exon + 1, y = max(cnvs_within_exon$exon_overlap_length, na.rm = TRUE) / 2,
           label = paste0("Minimum percentage overlap of recalled CNVs = ", min_recalled_within_exon),
           color = "#1f77b4", hjust = 0, size = 4)

dot_plot_cnvs_within_exon
