library(ggplot2)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(stringr)
library(scales)

# # Combine single and multi exon simulation data
# sim_table_combined <- bind_rows(sim_singleexon_table_combined, 
#                                 sim_multiexon_table_combined)

# Summarise p-values for simulation data
p_values_df <- sim_table_combined %>%
  group_by(cn, exon_type) %>%
  summarise(
    mean_p_value = mean(log_p_value, na.rm = TRUE),
    ci_lower = mean(log_p_value, na.rm = TRUE) - qt(0.975, df = sum(!is.na(log_p_value)) - 1) * sd(log_p_value, na.rm = TRUE) / sqrt(sum(!is.na(log_p_value))),
    ci_upper = mean(log_p_value, na.rm = TRUE) + qt(0.975, df = sum(!is.na(log_p_value)) - 1) * sd(log_p_value, na.rm = TRUE) / sqrt(sum(!is.na(log_p_value)))
  )

# Ensure cn includes all CN values as factors
p_values_df <- p_values_df %>%
  mutate(cn = factor(cn, levels = c(0, 1, 2, 3, 4)))

# Ensure Estimated.CN.of.CNV includes all CN values as factors and handle outliers
sample_cnv_calls_combined <- sample_cnv_calls_combined %>%
  mutate(cn = factor(Estimated.CN.of.CNV, levels = c(0, 1, 2, 3, 4)))

# Create a new table with each row representing a sample with call
# and create a new column with number of CNVs
samples_with_calls <- sample_cnv_calls_combined %>%
  group_by(sample, sample_id_trimmed, run, outlier_sample, Percent.Difference) %>%
  summarise(number_of_cnvs = n()  # Number of CNVs found in the sample
  ) %>%
  ungroup()

# Note : these are samples that passed sequencing QC
all_samples_combined <- all_samples_combined %>%
  left_join(samples_with_calls, by = c("sample", 
                                       "sample_id_trimmed", 
                                       "run", "outlier_sample", 
                                       "Percent.Difference")) %>%
  mutate(number_of_cnvs = coalesce(number_of_cnvs, 0))


# Calculate total samples per run and mean number across runs
samples_summary <- all_samples_combined %>%
  group_by(run) %>%
  summarise(num_samples = n_distinct(sample)) %>%
  ungroup()

# Print total number of samples per run. Note : these are samples that passed sequencing QC
print(samples_summary)
# write.csv(samples_summary, "num_samples_per_run.csv", row.names = FALSE)

# Print mean number of samples per run
cat("\nMean number of samples per run:", mean(samples_summary$num_samples), "\n")


# write.csv(all_samples_combined, "all_samples_combined.csv", 
#           row.names = FALSE)


# Distribution of CNV calls per sample (SHOWN IN THESIS!!)
# Helper function for statistical mode
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Filter data
plot_data <- subset(all_samples_combined, number_of_cnvs > 0)

# Pre-calculate statistics
cnv_mean   <- mean(plot_data$number_of_cnvs)
cnv_median <- median(plot_data$number_of_cnvs)
cnv_mode   <- get_mode(plot_data$number_of_cnvs)

# Plot
ggplot(plot_data, aes(x = number_of_cnvs)) +
  geom_histogram(bins = 40, fill = "blue", color = "black") +
  geom_vline(aes(xintercept = cnv_mean,   color = "Mean"),   linewidth = 1, linetype = "solid") +
  geom_vline(aes(xintercept = cnv_median, color = "Median"), linewidth = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = cnv_mode,   color = "Mode"),   linewidth = 1, linetype = "dotted") +
  annotate("text", x = cnv_mean,   y = Inf, vjust = 1.5, hjust = -0.1,
           label = paste("Mean =", round(cnv_mean, 1)),   color = "red") +
  annotate("text", x = cnv_median, y = Inf, vjust = 3.0, hjust = -0.1,
           label = paste("Median =", round(cnv_median, 1)), color = "orange") +
  annotate("text", x = cnv_mode,   y = Inf, vjust = 4.5, hjust = -0.1,
           label = paste("Mode =", cnv_mode),               color = "darkgreen") +
  scale_color_manual(
    name   = "Central Tendency",
    values = c(Mean = "red", Median = "orange", Mode = "darkgreen")
  ) +
  labs(
    title = "Distribution of VS-CNV calls per BAM for BAMs with >= 1 VS-CNV call",
    x     = "Number of VS-CNV calls",
    y     = "Number of BAMs"
  ) +
  theme_minimal()


# Distribution of CNV calls per sample for non-outliers
ggplot(all_samples_combined %>% filter(outlier_sample == 'NO'), aes(x = number_of_cnvs)) +
  geom_histogram(bins = 20, fill = "blue", color = "black") +
  scale_y_continuous(breaks = pretty_breaks()) +  # Ensure whole numbers on y-axis
  labs(title = "Distribution of CNV calls per sample for non-outliers",
       x = "Number of CNVs",
       y = "Count") +
  theme_minimal()

# Distribution of CNV calls per sample for outliers
ggplot(all_samples_combined %>% filter(outlier_sample == 'YES'), aes(x = number_of_cnvs)) +
  geom_histogram(bins = 20, fill = "blue", color = "black") +
  scale_y_continuous(breaks = pretty_breaks()) +  # Ensure whole numbers on y-axis
  labs(title = "Distribution of CNV calls per sample for outliers",
       x = "Number of CNVs",
       y = "Count") +
  theme_minimal()


# Correlating Percent Difference and Number of CNVs called
cor_coeff <- cor(all_samples_combined$number_of_cnvs, all_samples_combined$Percent.Difference, use = "complete.obs")

# Create the linear regression plot with the correlation coefficient
ggplot(all_samples_combined, aes(x = number_of_cnvs, y = Percent.Difference)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  # Plot the points
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Add the linear regression line
  labs(title = "Linear Regression of Percent Difference vs. Number of CNVs called",
       x = "Number of CNVs",
       y = "Percent Difference") +
  annotate("text", x = Inf, y = Inf, label = paste("r =", round(cor_coeff, 2)), 
           vjust = 2, hjust = 2, size = 5, color = "black") +  # Add the correlation coefficient as a label
  theme_minimal()


# Samples with multiple calls. Check if they are all duplications or deletions.
multiple_calls <- sample_cnv_calls_combined %>%
  group_by(sample_id_trimmed, Estimated.CN.of.CNV) %>%
  filter(n() > 1) %>%  # Keep only samples with more than one CNV call
  summarise(cnv_calls_count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = Estimated.CN.of.CNV, values_from = cnv_calls_count, values_fill = 0)

# Filter for samples with only CN = 1
cn_1_samples <- multiple_calls %>%
  filter(`1` > 0 & `3` == 0 & `4` == 0)

# Filter for samples with only CN = 3
cn_3_samples <- multiple_calls %>%
  filter(`3` > 0 & `1` == 0 & `4` == 0)

# Filter for samples with only CN = 4
cn_4_samples <- multiple_calls %>%
  filter(`4` > 0 & `1` == 0 & `3` == 0)


all_samples_has_calls <- all_samples_combined %>%
  filter(number_of_cnvs > 0)


# Plot percent diff vs. number_of_cnvs (flipped axes, no run information, y-axis spaced by 2)
correlation_plot <- ggplot(all_samples_has_calls, aes(x = Percent.Difference, y = number_of_cnvs)) +
  geom_point(size = 3, color = "blue") +
  labs(
    title = "Number of CNVs vs. Percent Difference",
    x = "Percent Difference",
    y = "Number of CNVs"
  ) +
  scale_y_continuous(breaks = seq(0, max(all_samples_has_calls$number_of_cnvs, na.rm = TRUE), by = 2)) +
  theme_minimal() +
  theme(
    panel.grid = element_line(color = "gray85"),
    legend.position = "none" # Remove legend since run is no longer included
  )

correlation_plot


