library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tibble)
library(scales)
library(plotly)
library(rlang)

### For top 10 VS-CNV calls comparison with runs on validation set ###

dir <- "/path/to/project/Aim_1/CNV_calling/cnvbenchmarker/CNVbenchmarkeR2_test/output"

# Create output folder
output_plot_dir <- file.path(dir, "zscore_density_plots")
if (!dir.exists(output_plot_dir)) dir.create(output_plot_dir)

# List of harmonised data frames and their quality score variables
harmonised_list <- list(
  "VS-CNV_perrun" = list(data = vscnv_reorder_perrun, quality_col = "Avg.Z.Score.of.CNV"),
  "VS-CNV_val" = list(data = vscnv_reorder_val, quality_col = "Avg.Z.Score.of.CNV"),
  "atlasCNV"   = list(data = atlascnv_reorder,     quality_col = "c_Score"),
  "clearCNV"   = list(data = clearcnv_reorder,     quality_col = "score"),
  "clinCNV"    = list(data = clincnv_reorder,      quality_col = "loglikelihood"),
  "cnvkit"     = list(data = cnvkit_reorder,       quality_col = "log2"),
  "cobalt"     = list(data = cobalt_reorder,       quality_col = "quality"),
  "codex2"     = list(data = codex2_reorder,       quality_col = "mBIC"),
  "decon"      = list(data = decon_reorder,        quality_col = "BF"),
  "exomedepth" = list(data = exomedepth_reorder,   quality_col = "BF"),
  "gatk"       = list(data = gatk_reorder,         quality_col = "gt_QS"),
  "VisCap"     = list(data = viscap_reorder,       quality_col = "quality"),
  "panelcn"     = list(data = panelcn_reorder,     quality_col = "None")
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
    df %>% select(Tool, Sample, Gene, CNV_Type, zscore, score_col_name, raw_score)
  )
}


# mark outlier
outlier_sample <- vscnv_reorder_perrun %>%
  filter(Sample == "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S108") %>%
  select(Sample, Gene, CNV_Type)

# Compute -log10(p-value), filter out the excluded sample, and select top 10 and 20 CNVs from 
# the VS-CNV per-run CNV list
top10_vscnv <- vscnv_reorder_perrun %>%
  filter(Sample != "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S108") %>%
  mutate(log10p = -log10(adjusted_p_value)) %>%
  arrange(desc(log10p)) %>%
  slice_head(n = 10) %>%
  select(Sample, Gene, CNV_Type)

top20_vscnv <- vscnv_reorder_perrun %>%
  filter(Sample != "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S108") %>%
  mutate(log10p = -log10(adjusted_p_value)) %>%
  arrange(desc(log10p)) %>%
  slice_head(n = 20) %>%
  select(Sample, Gene, CNV_Type)


# Make sure Tool is a factor with consistent ordering
zscore_combined$Tool <- factor(zscore_combined$Tool, levels = sort(unique(zscore_combined$Tool)))

# Join flags into zscore_combined
zscore_combined <- zscore_combined %>%
  left_join(outlier_sample %>% mutate(outlier = TRUE), by = c("Sample", "Gene", "CNV_Type")) %>%
  left_join(top10_vscnv %>% mutate(Top10 = TRUE), by = c("Sample", "Gene", "CNV_Type")) %>%
  left_join(top20_vscnv %>% mutate(Top20 = TRUE), by = c("Sample", "Gene", "CNV_Type")) %>%
  mutate(
    Highlight = case_when(
      !is.na(Top10) ~ "Top10",
      !is.na(Top20) ~ "Top11-20",
      !is.na(outlier) ~ "Outlier",
      TRUE ~ "Other"
    )
  )

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
      "CNV_Type: ", CNV_Type, "<br>",
      "Quality metric (", score_col_name, "): ", raw_score, "<br>",
      "Z-score: ", round(zscore, 2)
    )
  )

plotly_combined <- plot_ly(
  data = zscore_combined,
  x = ~x_jittered,
  y = ~zscore,
  type = "scatter",
  mode = "markers",
  color = ~Highlight,
  colors = c("Other" = "gray", "Outlier" = "red", "Top11-20" = "blue", "Top10" = "green"),  # target_sample shown in blue for distinction
  text = ~text,
  hoverinfo = "text",
  marker = list(size = 7, opacity = 0.6)
) %>%
  layout(
    title = "Z-scaled Quality Scores per Tool for Validation Set",
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

# # Save interactive HTML
# htmlwidgets::saveWidget(
#   plotly_combined,
#   file = file.path(output_plot_dir, "zscore_combined_tools_plotly_top1-20.html"),
#   selfcontained = TRUE
# )
