library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tibble)
library(scales)
library(plotly)
library(rlang)


### run the overlap_calls_from_tools_and_vscnv.R script to get the variables first ####
#### Plot z-scaled quality scores for each tool with a target sample ####
### For Run 31 only ###

dir <- "/path/to/project/Aim_1/CNV_calling/cnvbenchmarker/CNVbenchmarkeR2_test/output_run10"

# Target sample
target_sample <- "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S36"
target_gene <- "APOB"

# Create output folder
output_plot_dir <- file.path(dir, "zscore_density_plots")
if (!dir.exists(output_plot_dir)) dir.create(output_plot_dir)

# List of harmonised data frames and their quality score variables
harmonised_list <- list(
  "VS-CNV_perrun" = list(data = vscnv_reorder_perrun, quality_col = "adjusted_p_value", transform = function(x) -log10(x)),
  "atlasCNV"   = list(data = atlascnv_reorder,     quality_col = "c_Score"),
  "clearCNV"   = list(data = clearcnv_reorder,     quality_col = "score"),
  "clinCNV"    = list(data = clincnv_reorder,      quality_col = "loglikelihood"),
  "cnvkit"     = list(data = cnvkit_reorder,       quality_col = "log2"),
  "cobalt"     = list(data = cobalt_reorder,       quality_col = "quality"),
  "codex2"     = list(data = codex2_reorder,       quality_col = "mBIC"),
  "decon"      = list(data = decon_reorder,        quality_col = "BF"),
  "exomedepth" = list(data = exomedepth_reorder,   quality_col = "BF"),
  "gatk"       = list(data = gatk_reorder,         quality_col = "gt_QS"),
  "VisCap"     = list(data = viscap_reorder,       quality_col = "quality")
)
# Note: panelcn does not output any quality scores...

# Combined z-score data storage
zscore_combined <- data.frame()

# Loop through tools to collect z-scored data
for (tool_name in names(harmonised_list)) {
  tool_info <- harmonised_list[[tool_name]]
  df <- tool_info$data
  score_col <- tool_info$quality_col
  transform_fun <- if (!is.null(tool_info$transform)) tool_info$transform else function(x) x
  
  raw_scores <- transform_fun(df[[score_col]])
  valid_scores <- raw_scores[!is.na(raw_scores)]
  if (length(valid_scores) < 3) next
  
  z_scores <- scale(valid_scores)
  df$zscore <- NA
  df$zscore[!is.na(raw_scores)] <- as.numeric(z_scores)
  
  df$Tool <- tool_name
  df$score_col_name <- score_col  # store quality metric name
  df$raw_score <- df[[score_col]]  # store raw score itself
  
  zscore_combined <- bind_rows(
    zscore_combined,
    df %>% select(Tool, Sample, Gene, zscore, score_col_name, raw_score)
  )
}

# Flag the target point
zscore_combined <- zscore_combined %>%
  mutate(
    Highlight = ifelse(Sample == target_sample & Gene == target_gene, target_sample, "Other")
  )

# Plot all tools in one figure
combined_plot <- ggplot(zscore_combined, aes(x = Tool, y = zscore, color = Highlight)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  scale_color_manual(values = setNames(c("red", "grey40"), c(target_sample, "Other"))) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Z-scaled Quality Scores per Tool",
    x = "Tool",
    y = "Z-scaled Quality Score for Run 31",
    color = "CNV"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

# Save to file
ggsave(
  filename = file.path(output_plot_dir, "zscore_combined_tools_run31.png"),
  plot = combined_plot,
  width = 10,
  height = 5,
  dpi = 300,
  bg = "white"
)

### generate plotly ###

# Make sure Tool is a factor with consistent ordering
zscore_combined$Tool <- factor(zscore_combined$Tool, levels = sort(unique(zscore_combined$Tool)))

# Add jittered numeric x-coordinates for plotting
set.seed(123)  # for reproducibility
zscore_combined <- zscore_combined %>%
  mutate(
    Tool_num = as.numeric(factor(Tool)),  # ensure numeric conversion is consistent
    x_jittered = Tool_num + runif(n(), min = -0.2, max = 0.2),
    text = paste0(
      "Tool: ", Tool, "<br>",
      "Sample: ", Sample, "<br>",
      "Gene: ", Gene, "<br>",
      "Quality metric (", score_col_name, "): ", raw_score, "<br>",
      "Z-score: ", round(zscore, 2)
    )
  )

# Create interactive jittered scatter plot
plotly_combined <- plot_ly(
  data = zscore_combined,
  x = ~x_jittered,
  y = ~zscore,
  type = "scatter",
  mode = "markers",
  color = ~Highlight,
  colors = setNames(c("gray", "red"), c("Other", target_sample)),
  text = ~text,
  hoverinfo = "text",
  marker = list(size = 7, opacity = 0.6)
) %>%
  layout(
    title = "Z-scaled Quality Scores per Tool for Run 31",
    xaxis = list(
      title = "Tool",
      tickmode = "array",
      tickvals = 1:length(levels(zscore_combined$Tool)),
      ticktext = levels(zscore_combined$Tool),
      tickangle = -45
    ),
    yaxis = list(title = "Z-scaled Quality Score"),
    legend = list(title = list(text = "Sample"))
  )

# Save interactive HTML
htmlwidgets::saveWidget(
  plotly_combined,
  file = file.path(output_plot_dir, "zscore_combined_tools_plotly_run31.html"),
  selfcontained = TRUE
)


