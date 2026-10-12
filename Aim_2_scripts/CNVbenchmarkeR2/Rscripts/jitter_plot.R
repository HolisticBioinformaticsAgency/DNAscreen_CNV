library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(tibble)
library(scales)
library(plotly)
library(rlang)
library(IRanges)
library(htmlwidgets)
library(GenomicRanges)
library(patchwork)

# ============================================================
# Helper functions
# ============================================================

wilson_ci <- function(x, n, conf.level = 0.95) {
  if (n == 0) {
    return(c(
      estimate = NA_real_,
      lower = NA_real_,
      upper = NA_real_
    ))
  }
  
  z <- qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  denominator <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denominator
  margin <- (z / denominator) *
    sqrt((p * (1 - p) / n) + z^2 / (4 * n^2))
  
  c(
    estimate = p,
    lower = max(0, centre - margin),
    upper = min(1, centre + margin)
  )
}

# ============================================================
# Construct combined z-score data
# ============================================================

zscore_combined <- data.frame()

for (tool_name in names(harmonised_list_with_nonvscnv)) {
  tool_info <- harmonised_list_with_nonvscnv[[tool_name]]
  df <- tool_info$data
  
  if (tool_name == "panelcn") {
    df$Tool <- tool_name
    df$raw_score <- NA_real_
    df$score_col_name <- "None"
    df$zscore <- 0
  } else {
    score_col <- tool_info$quality_col
    raw_scores <- as.numeric(df[[score_col]])
    valid_scores <- abs(raw_scores[!is.na(raw_scores)])
    
    if (length(valid_scores) < 3) {
      next
    }
    
    df$zscore <- NA_real_
    df$zscore[!is.na(raw_scores)] <- as.numeric(scale(valid_scores))
    df$Tool <- tool_name
    df$score_col_name <- score_col
    df$raw_score <- raw_scores
  }
  
  zscore_combined <- bind_rows(
    zscore_combined,
    df %>%
      select(
        Tool, Sample, Gene, Chromosome, Start, End, CNV_Type,
        zscore, score_col_name, raw_score, acgh_status, non_vscnv
      )
  )
}

# ============================================================
# Non-VS-CNV presentation plot
# ============================================================

zscore_combined_no_vscnv <- zscore_combined %>%
  filter(!acgh_status %in% c("tbc", "inconclusive")) %>%
  filter(
    Tool != "VS-CNV_val",
    Gene != "ATM"
  ) %>%
  mutate(
    Tool = if_else(Tool == "VS-CNV_perrun", "VS-CNV", Tool),
    CNV_color = case_when(
      acgh_status == "concordant" ~ "Concordant with aCGH",
      acgh_status == "discordant" ~ "Discordant with aCGH",
      TRUE ~ NA_character_
    ),
    CNV_shape = if_else(
      non_vscnv == "yes",
      "Not called by VS-CNV (preliminary CNV calling)",
      "Called by VS-CNV (preliminary CNV calling)"
    )
  )

tool_order_no_vscnv <- c("VS-CNV", other_tools)

zscore_combined_no_vscnv$Tool <- factor(
  zscore_combined_no_vscnv$Tool,
  levels = tool_order_no_vscnv
)

zscore_combined_no_vscnv <- zscore_combined_no_vscnv %>%
  filter(Tool != "VS-CNV")

combined_plot_no_vscnv <- ggplot(
  zscore_combined_no_vscnv,
  aes(x = Tool, y = zscore)
) +
  geom_jitter(
    aes(color = CNV_color, shape = CNV_shape),
    width = 0.2,
    alpha = 0.6,
    size = 1.5
  ) +
  scale_color_manual(
    values = c(
      "Concordant with aCGH" = "green",
      "Discordant with aCGH" = "red"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Not called by VS-CNV (preliminary CNV calling)" = 16,
      "Called by VS-CNV (preliminary CNV calling)" = 16
    )
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Z-scaled quality scores per tool",
    x = "Tool",
    y = "Z-scaled quality score",
    color = "aCGH status",
    shape = "Preliminary CNV calling"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

combined_plot_no_vscnv

# ============================================================
# Interactive non-VS-CNV plot
# ============================================================

combined_plotly_no_vscnv <- ggplot(
  zscore_combined_no_vscnv,
  aes(
    x = Tool,
    y = zscore,
    color = CNV_color,
    shape = CNV_shape,
    text = paste(
      "Sample:", Sample,
      "<br>Gene:", Gene,
      "<br>Coordinates:", paste0("chr", Chromosome, ":", Start, "-", End),
      "<br>Raw score:", round(raw_score, 3),
      "<br>Z-score:", round(zscore, 2),
      "<br>Status:", CNV_color
    )
  )
) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c(
      "Concordant with aCGH" = "green",
      "Discordant with aCGH" = "red"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Not called by VS-CNV (preliminary CNV calling)" = 17,
      "Called by VS-CNV (preliminary CNV calling)" = 16
    )
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Z-scaled quality scores per tool",
    x = "Tool",
    y = "Z-scaled quality score",
    color = "aCGH status",
    shape = "Preliminary CNV calling"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

interactive_plot_no_vscnv <- ggplotly(
  combined_plotly_no_vscnv,
  tooltip = "text"
)

# interactive_plot_no_vscnv

# ============================================================
# Thesis plot including VS-CNV
# ============================================================

zscore_combined_with_vscnv <- zscore_combined %>%
  filter(Gene != "ATM") %>%
  filter(Tool != "VS-CNV_perrun") %>%
  mutate(
    Tool = if_else(Tool == "VS-CNV_val", "VS-CNV", Tool)
  ) %>%
  filter(!acgh_status %in% c("inconclusive", "tbc")) %>%
  mutate(
    CNV_color = case_when(
      acgh_status == "concordant" ~ "True-positive Call",
      acgh_status == "discordant" ~ "False-positive Call",
      TRUE ~ NA_character_
    ),
    CNV_shape = if_else(
      non_vscnv == "yes",
      "Not called by VS-CNV (preliminary CNV calling)",
      "Called by VS-CNV (preliminary CNV calling)"
    )
  )

tool_order_with_vscnv <- c(other_tools, "VS-CNV")

zscore_combined_with_vscnv$Tool <- factor(
  zscore_combined_with_vscnv$Tool,
  levels = tool_order_with_vscnv
)

tool_labels <- c(
  "VS-CNV" = "VS-CNV",
  "atlasCNV" = "Atlas-CNV",
  "clearCNV" = "ClearCNV",
  "clinCNV" = "ClinCNV",
  "cnvkit" = "CNVkit",
  "cobalt" = "Cobalt",
  "codex2" = "CODEX2",
  "decon" = "DECoN",
  "exomedepth" = "ExomeDepth",
  "gatk" = "GATK-gCNV",
  "VisCap" = "VisCap",
  "panelcn" = "panelcn.MOPS"
)

# ============================================================
# Read validated truth-set regions
# ============================================================

validated_res <- read.csv(
  "/path/to/project/Aim_1/Validations/acgh_validation_status_files/acgh_round_one_res_for_candidate_cnvs_dan.csv"
)

validated_res$Man..cur..Chr <- gsub(
  "^chr",
  "",
  validated_res$Man..cur..Chr
)

validated_res$Var.Seq.Chr <- gsub(
  "^chr",
  "",
  validated_res$Var.Seq.Chr
)

validated_res[] <- lapply(validated_res, function(x) {
  if (is.character(x)) {
    x[x == ""] <- NA
  }
  x
})

tp_truth <- validated_res %>%
  filter(Verified_20250730 == "YES") %>%
  filter(
    !is.na(Man..cur..Chr),
    !is.na(Manual.curation.Start),
    !is.na(Manual.curation.Stop)
  ) %>%
  mutate(
    tp_region_id = paste0("TPregion_", row_number()),
    composite_chr = paste0(sample, "_", Man..cur..Chr)
  )

tp_gr <- GRanges(
  seqnames = tp_truth$composite_chr,
  ranges = IRanges(
    start = as.numeric(tp_truth$Manual.curation.Start),
    end = as.numeric(tp_truth$Manual.curation.Stop)
  ),
  tp_region_id = tp_truth$tp_region_id
)

# ============================================================
# Match caller outputs to truth-set regions
# ============================================================

tp_calls_df <- zscore_combined_with_vscnv %>%
  filter(CNV_color == "True-positive Call") %>%
  filter(
    !is.na(Chromosome),
    !is.na(Start),
    !is.na(End)
  ) %>%
  mutate(
    Chromosome = gsub("^chr", "", Chromosome),
    composite_chr = paste0(Sample, "_", Chromosome)
  )

calls_gr <- GRanges(
  seqnames = tp_calls_df$composite_chr,
  ranges = IRanges(
    start = as.numeric(tp_calls_df$Start),
    end = as.numeric(tp_calls_df$End)
  )
)

hits <- findOverlaps(calls_gr, tp_gr, minoverlap = 1L)

tp_calls_df$tp_region_id <- NA_character_
tp_calls_df$tp_region_id[queryHits(hits)] <-
  tp_gr$tp_region_id[subjectHits(hits)]

# ============================================================
# Calculate TP, FP, recall, precision, and Wilson intervals
# ============================================================

n_truth_positive <- 10

tp_counts <- tp_calls_df %>%
  group_by(Tool) %>%
  summarise(
    n_concordant = n_distinct(tp_region_id, na.rm = TRUE),
    .groups = "drop"
  )

fp_counts <- zscore_combined_with_vscnv %>%
  group_by(Tool) %>%
  summarise(
    n_discordant = sum(acgh_status == "discordant", na.rm = TRUE),
    .groups = "drop"
  )

all_tools <- tibble(
  Tool = factor(
    levels(zscore_combined_with_vscnv$Tool),
    levels = levels(zscore_combined_with_vscnv$Tool)
  )
)

perf_counts <- all_tools %>%
  left_join(tp_counts, by = "Tool") %>%
  left_join(fp_counts, by = "Tool") %>%
  mutate(
    n_concordant = replace_na(n_concordant, 0),
    n_discordant = replace_na(n_discordant, 0),
    n_total_calls = n_concordant + n_discordant
  ) %>%
  rowwise() %>%
  mutate(
    recall_est = wilson_ci(n_concordant, n_truth_positive)["estimate"],
    recall_lo = wilson_ci(n_concordant, n_truth_positive)["lower"],
    recall_hi = wilson_ci(n_concordant, n_truth_positive)["upper"],
    precision_est = wilson_ci(n_concordant, n_total_calls)["estimate"],
    precision_lo = wilson_ci(n_concordant, n_total_calls)["lower"],
    precision_hi = wilson_ci(n_concordant, n_total_calls)["upper"]
  ) %>%
  ungroup() %>%
  mutate(
    label_tp = paste0("N = ", n_concordant),
    label_fp = paste0("N = ", n_discordant),
    value_recall = sprintf(
      "%.0f%% (%.0f–%.0f%%)",
      recall_est * 100,
      recall_lo * 100,
      recall_hi * 100
    ),
    value_precision = if_else(
      n_total_calls > 0,
      sprintf(
        "%.0f%% (%.0f–%.0f%%)",
        precision_est * 100,
        precision_lo * 100,
        precision_hi * 100
      ),
      "NA"
    )
  )

# ============================================================
# Layout values
# ============================================================

# Visible jitter plot ends at Z = -2.5.
y_panel_min <- -2.5

# Summary rows are intentionally below the visible panel.
y_tp <- -2.95
y_fp <- -3.30
y_recall <- -3.75
y_precision <- -4.20

# Caller labels sit below the summary rows.
y_tool_name <- -5.45

perf_counts <- perf_counts %>%
  mutate(
    y_tp = y_tp,
    y_fp = y_fp,
    y_recall = y_recall,
    y_precision = y_precision
  )

tool_name_data <- tibble(
  Tool = factor(
    levels(zscore_combined_with_vscnv$Tool),
    levels = levels(zscore_combined_with_vscnv$Tool)
  ),
  tool_label = tool_labels[levels(zscore_combined_with_vscnv$Tool)],
  y_tool_name = y_tool_name
)

# ============================================================
# TP/FP labels centred around Z = 8
# ============================================================

y_tp_label <- 8.3
y_fp_label <- 7.7

perf_counts <- perf_counts %>%
  mutate(
    label_tp = paste0("TP = ", n_concordant),
    label_fp = paste0("FP = ", n_discordant)
  )


# ============================================================
# Reshape performance metrics for bar panel
# ============================================================

perf_long <- perf_counts %>%
  pivot_longer(
    cols = c(
      recall_est, recall_lo, recall_hi,
      precision_est, precision_lo, precision_hi
    ),
    names_to = c("metric", ".value"),
    names_pattern = "(recall|precision)_(est|lo|hi)"
  ) %>%
  mutate(
    metric = factor(
      recode(metric, recall = "Recall", precision = "Precision"),
      levels = c("Recall", "Precision")
    )
  )

# Tools with zero calls have no precision estimate
na_precision <- perf_long %>%
  filter(metric == "Precision", is.na(est))

# ============================================================
# Vertical dividers between callers
# ============================================================

n_tools <- nlevels(zscore_combined_with_vscnv$Tool)
caller_boundaries <- seq(1.5, n_tools - 0.5, by = 1)

# ============================================================
# Panel A: jitter plot (summary rows removed)
# ============================================================

set.seed(123)

jitter_panel <- ggplot(
  zscore_combined_with_vscnv,
  aes(x = Tool, y = zscore)
) +
  geom_jitter(
    aes(color = CNV_color),
    width = 0.2,
    alpha = 0.6,
    size = 1.5
  ) +
  geom_vline(
    xintercept = caller_boundaries,
    linetype = "dotted",
    color = "grey60",
    linewidth = 0.4
  ) +
  # TP and FP call counts above each caller
  geom_text(
    data = perf_counts,
    aes(x = Tool, y = y_tp_label, label = label_tp),
    color = "#2196F3",
    size = 3.5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = perf_counts,
    aes(x = Tool, y = y_fp_label, label = label_fp),
    color = "#F44336",
    size = 3.5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c(
      "True-positive Call" = "#2196F3",
      "False-positive Call" = "#F44336"
    ),
    breaks = c("True-positive Call", "False-positive Call"),
    na.translate = FALSE
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Jitter plot of Z-scaled caller-specific confidence metrics and bar plot of recall and precision",
    x = NULL,
    y = "Z-scaled confidence metric",
    color = "aCGH status:"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 10, r = 15, b = 0, l = 15)
  ) +
  coord_cartesian(
    ylim = c(-2.5, max(9, max(zscore_combined_with_vscnv$zscore, na.rm = TRUE) + 0.5))
  )

# ============================================================
# Panel B: recall and precision bars with Wilson 95% CI
# ============================================================

metric_colours <- c("Recall" = "#D46600", "Precision" = "#4A90E2")

bar_panel <- ggplot(
  perf_long %>% filter(!is.na(est)),
  aes(x = Tool, y = est, fill = metric)
) +
  geom_vline(
    xintercept = caller_boundaries,
    linetype = "dotted",
    color = "grey60",
    linewidth = 0.4
  ) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.5
  ) +
  geom_text(
    aes(y = hi, label = percent(est, accuracy = 1)),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3
  ) +
  # Flag tools with no calls (precision undefined)
  geom_text(
    data = na_precision,
    aes(x = Tool, y = 0.02, label = "NA"),
    inherit.aes = FALSE,
    size = 3,
    color = "grey40"
  ) +
  scale_fill_manual(values = metric_colours) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.08),
    expand = c(0, 0)
  ) +
  scale_x_discrete(labels = tool_labels) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Caller",
    y = "Percentage (Recall/Precision)",
    fill = "Metric:"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 0, r = 15, b = 10, l = 15)
  )

# ============================================================
# Stack panels and save
# ============================================================

combined_figure <- jitter_panel / bar_panel +
  plot_layout(heights = c(2.2, 1), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

combined_figure

# ggsave(
#   filename = file.path(
#     "/path/to/project/Thesis/thesis_aim_2_chapter/important_figures/zscore_combined_tools_acgh_thesis_with_vscnv_v3.png"
#   ),
#   plot = combined_figure,
#   width = 15,
#   height = 12,
#   dpi = 300,
#   bg = "white"
# )


# ============================================================
# Interactive plot with VS-CNV
# ============================================================

combined_plotly_with_vscnv <- ggplot(
  zscore_combined_with_vscnv,
  aes(
    x = Tool,
    y = zscore,
    color = CNV_color,
    shape = CNV_shape,
    text = paste(
      "Sample:", Sample,
      "<br>Gene:", Gene,
      "<br>Coordinates:", paste0("chr", Chromosome, ":", Start, "-", End),
      "<br>Raw score:", round(raw_score, 3),
      "<br>Z-score:", round(zscore, 2),
      "<br>Status:", CNV_color
    )
  )
) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c(
      "True-positive Call" = "#2196F3",
      "False-positive Call" = "#F44336"
    ),
    na.translate = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "Not called by VS-CNV (preliminary CNV calling)" = 17,
      "Called by VS-CNV (preliminary CNV calling)" = 16
    )
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Z-scaled caller-specific confidence metrics",
    x = "Caller",
    y = "Z-scaled quality score",
    color = "aCGH status",
    shape = "Preliminary CNV calling"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

interactive_plot_with_vscnv <- ggplotly(
  combined_plotly_with_vscnv,
  tooltip = "text"
)

# interactive_plot_with_vscnv