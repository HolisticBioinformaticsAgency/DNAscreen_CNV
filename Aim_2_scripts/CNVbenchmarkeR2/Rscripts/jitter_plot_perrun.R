library(dplyr)
library(ggplot2)
library(plotly)
library(tibble)
library(htmlwidgets)

dir <- "/path/to/project/Aim_1/CNV_calling/cnvbenchmarker/CNVbenchmarkeR2_test/output"


# Create output folder
output_plot_dir <- file.path(dir, "zscore_density_plots")
if (!dir.exists(output_plot_dir)) dir.create(output_plot_dir)

# --- Tools list (3 tools you want plotted) ---
tools_list <- list(
  "DECoN"   = decon_with_acgh,
  "GATK"    = gatk_with_acgh,
  "ClearCNV"= clearcnv_with_acgh
)

# Combined z-score data storage
zscore_combined <- tibble()

for (tool_name in names(tools_list)) {
  df <- tools_list[[tool_name]] %>% as.data.frame(stringsAsFactors = FALSE)
  
  # Ensure some standard columns exist and preserve originals
  # Use 'Quality' as the per-tool quality column (created earlier in your standardise_*)
  if (!"Quality" %in% names(df)) {
    # try to fall back on other column names if present
    if ("quality_col" %in% names(df) && !is.null(df[[df$quality_col[1]]])) {
      df$Quality <- df[[ df$quality_col[1] ]]
    } else if ("score" %in% names(df)) {
      df$Quality <- df$score
    } else {
      df$Quality <- NA_real_
    }
  }
  
  # Harmonise CNV type column name
  if (!"CNV_Type" %in% names(df)) {
    if ("CNV_type" %in% names(df)) df$CNV_Type <- df$CNV_type
    else if ("CNV.Type" %in% names(df)) df$CNV_Type <- df[["CNV.Type"]]
    else if ("cn" %in% names(df)) df$CNV_Type <- df$cn
    else df$CNV_Type <- NA_character_
  }
  
  # Ensure non_vscnv exists (may be missing)
  if (!"non_vscnv" %in% names(df)) df$non_vscnv <- NA_character_
  if (!"acgh_status" %in% names(df)) df$acgh_status <- NA_character_
  
  # Ensure Sample, Gene, Chromosome, Start, End columns exist (if not, create placeholders)
  for (col in c("Sample","Gene","Chromosome","Start","End","sample_id_trimmed")) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  
  # compute zscore per tool based on abs(Quality) for non-NA values
  raw_scores <- suppressWarnings(as.numeric(df$Quality))
  valid_idx <- which(!is.na(raw_scores))
  
  if (length(valid_idx) >= 3) {
    valid_scores <- abs(raw_scores[valid_idx])
    z_scores <- as.numeric(scale(valid_scores))
    df$zscore <- NA_real_
    df$zscore[valid_idx] <- z_scores
  } else {
    # not enough values to z-score — keep zscore as NA for all rows
    df$zscore <- NA_real_
  }
  
  # attach metadata columns similar to your original pipeline
  df$Tool <- tool_name
  df$score_col_name <- "Quality"
  df$raw_score <- raw_scores
  
  # Keep the columns your plotting code expects (add missing as NA)
  keep_cols <- c("Tool", "Sample", "run", "Gene", "Chromosome", "Start", "End",
                 "CNV_Type", "zscore", "score_col_name", "raw_score",
                 "acgh_status", "non_vscnv", "sample_id_trimmed")
  for (k in keep_cols) if (!k %in% names(df)) df[[k]] <- NA
  
  zscore_combined <- bind_rows(zscore_combined, df[keep_cols])
}

# Optional: remove ATM genes as in your reference
zscore_combined <- zscore_combined %>% filter(Gene != "ATM")

# Assign CNV color/shape categories (as in your reference)
sample_with_unconfirmed_cnv = c("DNS-XTHS-xxxx-xxx-DNSxxxxxx_S225",
                                "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S182",
                                "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S317",
                                "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S238",
                                "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S115",
                                "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S3",
                                "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S224",
                                "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S205")

zscore_combined <- zscore_combined %>%
  mutate(
    CNV_color = case_when(
      acgh_status == "concordant" ~ "Concordant with aCGH",
      acgh_status == "discordant" ~ "Discordant with aCGH",
      Sample %in% sample_with_unconfirmed_cnv ~ "From a Sample with tbc CNV",
      acgh_status == "tbc" ~ "From a Sample not on array yet",
      TRUE ~ "Unknown"
    ),
    CNV_shape = if_else(non_vscnv == "yes", "Non-VS-CNV", "VS-CNV")
  )

# Set tool order to control x-axis ordering
tool_order <- c("DECoN","GATK","ClearCNV")
zscore_combined$Tool <- factor(zscore_combined$Tool, levels = tool_order)

# Drop rows without any zscore if you want the plot to show only scored calls,
# otherwise comment out the next line to include NA zscores.
# zscore_combined <- zscore_combined %>% filter(!is.na(zscore))

# --- ggplot jitter (static) ---
jitter_plot <- ggplot(zscore_combined,
                      aes(x = Tool, y = zscore,
                          color = CNV_color, shape = CNV_shape,
                          text = paste("Sample:", Sample,
                                       "<br>Run:", run,
                                       "<br>Gene:", Gene,
                                       "<br>Coordinates:", paste0("chr", Chromosome, ":", Start, "-", End),
                                       "<br>Raw score:", round(raw_score, 3),
                                       "<br>Z-score:", round(zscore, 2),
                                       "<br>Status:", CNV_color))) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Concordant with aCGH" = "green",
    "Discordant with aCGH" = "red",
    "From a Sample with tbc CNV" = "pink",
    "From a Sample not on array yet" = "#666666",
    "Unknown" = "black"
  )) +
  scale_shape_manual(values = c("Non-VS-CNV" = 17, "VS-CNV" = 16)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Z-scaled Quality Scores per Tool",
    x = "Tool",
    y = "Z-scaled Quality Score",
    color = "aCGH status",
    shape = "CNV origin"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

# Show static
print(jitter_plot)

# --- Interactive version (plotly) ---
interactive_plot <- ggplotly(jitter_plot, tooltip = "text")
interactive_plot

# saveWidget(interactive_plot, file = file.path(output_plot_dir, "zscore_combined_tools_round_one_acgh_perrun.html"), selfcontained = TRUE)
