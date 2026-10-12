library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(forcats)
library(gridExtra)
library(grid)

########## Make sure to run vs-cnv_accuracy_metrics.R first before running this script!


dir <- "/path/to/project/Aim_1/CNV_calling/cnvbenchmarker/CNVbenchmarkeR2_test/output"

results_table <- read.csv(paste0(dir, "/summary/results_val_round_one.tsv"), sep = "\t")

results_table <- results_table %>%
  filter(algorithm != "VS-CNV_val",
         algorithm != "VS-CNV_perrun") %>%
  mutate(algorithm = recode(algorithm,
                              "VS-CNV_perrun" = "VS-CNV"))

####### TP and FN barplots at ROI and Gene Level #######

roi_lvl <- results_table %>% filter(level == "roi")

# Correct: Calculate max_TP as TP + FN for each tool (i.e., the total actual positives)
roi_tp_df <- roi_lvl %>%
  select(algorithm, TP, FN) %>%
  mutate(
    max_TP = TP + FN,
    diff = max_TP - TP,
    algorithm = factor(algorithm, levels = roi_lvl %>% arrange(desc(TP)) %>% pull(algorithm))
  ) %>%
  pivot_longer(cols = c(TP, diff), names_to = "segment", values_to = "value")

# Now the same fix for gene-level
gene_lvl <- results_table %>% filter(level == "gene")

gene_tp_df <- gene_lvl %>%
  select(algorithm, TP, FN) %>%
  mutate(
    max_TP = TP + FN,
    diff = max_TP - TP,
    algorithm = factor(algorithm, levels = gene_lvl %>% arrange(desc(TP)) %>% pull(algorithm))
  ) %>%
  pivot_longer(cols = c(TP, diff), names_to = "segment", values_to = "value")


# Plot with legend for extracting
p_roi_legend <- ggplot(roi_tp_df, aes(x = algorithm, y = value, fill = segment)) +
  geom_col() +
  scale_fill_manual(
    values = c("TP" = "#1b9e77", "diff" = "#d95f02"),
    labels = c("Number of aCGH-validated calls that are missed", "Number of aCGH-validated calls that are detected"),
    guide = guide_legend(title = NULL)
  ) +
  labs(
    x = "Tool",
    y = "Concordant calls",
    title = paste0("ROI Level ", "(N = ", roi_tp_df$max_TP[1], ")")
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

# Extract legend as a separate plot
legend <- cowplot::get_legend(p_roi_legend)

# Remove legend from original plots
p_roi <- p_roi_legend + theme(legend.position = "none")

p_gene <- ggplot(gene_tp_df, aes(x = algorithm, y = value, fill = segment)) +
  geom_col() +
  scale_fill_manual(
    values = c("TP" = "#1b9e77", "diff" = "#d95f02"),
    labels = c("Number of concordant calls that are missed", "Number of concordant calls that are detected"),
    guide = guide_legend(title = NULL)
  ) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 10), labels = scales::number_format(accuracy = 1)) +
  labs(
    x = "Tool",
    y = NULL,
    title = paste0("Maximum Concordant Calls ", "= ", gene_tp_df$max_TP[1])
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16), 
        legend.position = "none",
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16))


# Combine legend and plots: legend on top centered, plots side by side below
title <- textGrob("Number of concordant (aCGH-validated) \nCNV calls that are detected by each tool",
                  gp = gpar(fontsize = 16, fontface = "bold"))

p_combined <- grid.arrange(
  title,
  legend,
  gridExtra::arrangeGrob(p_gene, ncol = 1),
  nrow = 3,
  heights = c(0.1, 0.1, 1)
)

# Alternatively, if you want to stick with patchwork, you can do:
# p_combined <- wrap_elements(full = legend) / (p_roi | p_gene) + plot_layout(heights = c(0.1, 1))

# Print combined plot
print(p_combined)



# Save to file
# ggsave(
#   filename = file.path(dir, "summary/TP_FN_barplot_gene_lvl.png"),
#   plot = p_combined,
#   width = 8,
#   height = 5,
#   dpi = 300,
#   bg = "white"
# )




######## FP barplots for each tool at ROI and Gene Level ######

# Prepare ROI data for Discordant calls (False positives here)
roi_lvl <- results_table %>% filter(level == "roi")

roi_fp_df <- roi_lvl %>%
  select(algorithm, FP) %>%
  mutate(
    algorithm = factor(algorithm, levels = roi_lvl %>% arrange(FP) %>% pull(algorithm)),
    max_FP = max(FP, na.rm = TRUE),
    diff = FP  # We only care about FP here; no diff needed but keep for compatibility
  ) %>%
  pivot_longer(cols = c(FP), names_to = "segment", values_to = "value")

# Prepare Gene data for False Positives
gene_lvl <- results_table %>% filter(level == "gene")

gene_fp_df <- gene_lvl %>%
  select(algorithm, FP) %>%
  mutate(
    algorithm = factor(algorithm, levels = gene_lvl %>% arrange(FP) %>% pull(algorithm)),
    max_FP = max(FP, na.rm = TRUE),
    diff = FP
  ) %>%
  pivot_longer(cols = c(FP), names_to = "segment", values_to = "value")

# Plot for ROI level (FP)
p_roi_legend <- ggplot(roi_fp_df, aes(x = algorithm, y = value, fill = segment)) +
  geom_col() +
  scale_fill_manual(
    values = c("FP" = "red"),  # Pink/purple-ish for FPs
    labels = c("Number of discordant calls"),
    guide = guide_legend(title = NULL)
  ) +
  labs(
    x = "Tool",
    y = "Discordant calls",
    title = "ROI Level"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

legend <- cowplot::get_legend(p_roi_legend)

p_roi <- p_roi_legend + theme(legend.position = "none")

# Plot for Gene level (FP)
p_gene <- ggplot(gene_fp_df, aes(x = algorithm, y = value, fill = segment)) +
  geom_col() +
  scale_fill_manual(
    values = c("FP" = "red"),
    labels = c("Number of discordant calls"),
    guide = guide_legend(title = NULL)
  ) +
  labs(
    x = "Tool",
    y = NULL
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16), 
        legend.position = "none",
        axis.title.x = element_text(size = 16))

# Title grob
title <- textGrob("Number of discordant (NOT aCGH-validated) calls\n that are detected by each tool",
                  gp = gpar(fontsize = 16, fontface = "bold"))

# Combine legend and plots side by side with title on top
p_combined_fp <- grid.arrange(
  title,
  legend,
  gridExtra::arrangeGrob(p_gene, ncol = 1),
  nrow = 3,
  heights = c(0.1, 0.1, 1)
)

print(p_combined_fp)

# Save to file
# ggsave(
#   filename = file.path(dir, "summary/FP_barplot_gene_lvl.png"),
#   plot = p_combined_fp,
#   width = 8,
#   height = 5,
#   dpi = 300,
#   bg = "white"
# )
