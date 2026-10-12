library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tibble)
library(scales)
library(plotly)
library(rlang)
library(IRanges)
library(plotly)
library(htmlwidgets)

#### For HGSA poster!!! ###

dir <- "/path/to/project/Aim_1/CNV_calling/cnvbenchmarker/CNVbenchmarkeR2_test/output"


# Create output folder
output_plot_dir <- file.path(dir, "zscore_density_plots")
if (!dir.exists(output_plot_dir)) dir.create(output_plot_dir)


# # Create zeros for panelcn since it doesnt have quality scores..
# harmonised_list_with_acgh$panelcn$quality_col = "ones"
# harmonised_list_with_acgh$panelcn$data <- harmonised_list_with_acgh$panelcn$data %>%
#   mutate(ones = 1)

# Combined z-score data storage
zscore_combined <- data.frame()


# Loop through tools to collect z-scored data
for (tool_name in names(harmonised_list_with_acgh)) {
  tool_info <- harmonised_list_with_acgh[[tool_name]]
  df <- tool_info$data
  if (tool_name == "panelcn"){
    df$Tool <- tool_name
    df$raw_score <- NA
    df$score_col_name <- "None"
    df$zscore <- 0
  }
  else {
    score_col <- tool_info$quality_col
    # transform_fun <- if (!is.null(tool_info$transform)) tool_info$transform else function(x) x
    raw_scores <- as.numeric(df[[score_col]])
    valid_scores <- raw_scores[!is.na(raw_scores)]
    valid_scores <- abs(valid_scores)
    if (length(valid_scores) < 3) next
    
    z_scores <- scale(valid_scores)
    df$zscore <- NA
    df$zscore[!is.na(raw_scores)] <- as.numeric(z_scores)
    
    df$Tool <- tool_name
    df$score_col_name <- score_col  # store quality metric name
    df$raw_score <- df[[score_col]]  # store raw score itself
  }
  
  zscore_combined <- bind_rows(
    zscore_combined,
    df %>% select(Tool, Sample, Gene, Chromosome, Start, End, 
                  CNV_Type, zscore, score_col_name, raw_score, acgh_status)
  )
}



# Make sure Tool is a factor with consistent ordering
# zscore_combined$Tool <- factor(zscore_combined$Tool, levels = sort(unique(zscore_combined$Tool)))

### For presentations... we won't show the CNVs with tbc acgh status
zscore_combined <- zscore_combined %>%
  filter(acgh_status != "tbc")

# Set highlight for plotting and also filter out those we have yet
# to confirm aCGH validation status
zscore_combined <- zscore_combined %>%
  filter(Tool != "VS-CNV_val") %>%
  mutate(Tool = if_else(Tool == "VS-CNV_perrun", "VS-CNV", Tool)) %>%
  mutate(
    Highlight = case_when(
      acgh_status == "concordant" & Gene == "BRCA1" &
        Start <= 43045827 & End >= 43045653 &
        Sample == "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S96" ~ 
        "BRCA1 deletion concordant with aCGH",
      
      # acgh_status == "concordant" & Gene == "PCSK9" &
      #   Start <= 55040069 & End >= 55039813 &
      #   Sample == "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S210" ~ 
      #   "PCSK9 duplication concordant with aCGH",
      
      acgh_status == "discordant" & Gene == "APOB" &
        Start <= 21043970 & End >= 21037933 &
        Sample == "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S220" ~ 
        "APOB duplication discordant with aCGH",
      
      # acgh_status == "discordant" & Gene == "PCSK9" &
      #   Start <= 55063609 & End >= 55039813 &
      #   Sample == "DNS-XTHS-xxxx-xxx-DNSxxxxxx_S111" ~ 
      #   "PCSK9 duplication discordant with aCGH",
      acgh_status == "concordant" ~ "Concordant with aCGH",
      acgh_status == "discordant" ~ "Discordant with aCGH",
      acgh_status == "tbc" ~ "From a Sample with tbc aCGH status",
      TRUE ~ NA_character_
    )
  )

# 2️⃣ Add empty rows for Spacer1, aCGH, Spacer2
empty_cols <- tibble(
  Tool = c("Spacer1", "aCGH", "Spacer2"),
  Sample = NA_character_,
  Gene = NA_character_,
  Chromosome = NA_character_,
  Start = NA_real_,
  End = NA_real_,
  CNV_Type = NA_character_,
  zscore = NA_real_,
  score_col_name = NA_character_,
  raw_score = NA_real_,
  acgh_status = NA_character_,
  Highlight = NA_character_
)

zscore_combined <- bind_rows(zscore_combined, empty_cols)

# 3️⃣ Apply the consistent factor ordering
zscore_combined$Tool <- factor(zscore_combined$Tool, levels = tool_order)


# Filter for main violin categories (exclude tbc and others)
plot_data <- zscore_combined %>%
  filter(Highlight %in% c("Concordant with aCGH", "Discordant with aCGH") | is.na(Highlight))

# Special points data but mapped to main violin category for positioning
# special_points <- zscore_combined %>%
#   filter(Highlight %in% c("BRCA1 deletion concordant with aCGH",
#                           "APOB duplication discordant with aCGH")) %>%
#   mutate(Highlight_main = case_when(
#     Highlight == "BRCA1 deletion concordant with aCGH" ~ "Concordant with aCGH",
#     Highlight == "APOB duplication discordant with aCGH" ~ "Discordant with aCGH"
#   ))

special_points <- zscore_combined %>%
  filter(Highlight %in% c("BRCA1 deletion concordant with aCGH",
                          "APOB duplication discordant with aCGH")) %>%
  filter(!Tool %in% c("codex2", "exomedepth", "cnvkit")) %>%   # Exclude these tools here
  mutate(Highlight_main = case_when(
    Highlight == "BRCA1 deletion concordant with aCGH" ~ "Concordant with aCGH",
    Highlight == "APOB duplication discordant with aCGH" ~ "Discordant with aCGH"
  ))

# Factor levels
plot_data$Highlight <- factor(plot_data$Highlight,
                              levels = c("Concordant with aCGH", "Discordant with aCGH"))

special_points$Highlight_main <- factor(
  special_points$Highlight_main,
  levels = levels(plot_data$Highlight)
)

# Violin plot
violin_plot <- ggplot(plot_data, aes(x = Tool, y = zscore, fill = Highlight)) +
  geom_violin(width = 2, position = position_dodge(width = 1.0), alpha = 0.7, trim = FALSE) +
  # Overlay special points on correct violin position
  geom_point(
    data = special_points,
    aes(
      x = Tool,
      y = zscore,
      color = Highlight,
      shape = Highlight_main,
      group = Highlight_main,  # match to the violin grouping
      fill = Highlight_main    # add this so dodge matches violin fill groups
    ),
    position = position_dodge(width = 1.0),
    size = 3
  ) +
  scale_fill_manual(values = c("Concordant with aCGH" = "green",
                               "Discordant with aCGH" = "red")) +
  scale_color_manual(values = c("BRCA1 deletion concordant with aCGH" = "#7570B3",
                                "APOB duplication discordant with aCGH" = "#E6AB02")) +
  scale_shape_manual(values = c("Concordant with aCGH" = 16,
                                "Discordant with aCGH" = 16)) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Tool",
    y = "Z-scaled Quality Score",
    fill = "CNV Concordance",
    color = "Special CNV",
    shape = "CNV Concordance"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

violin_plot


ggsave(
  filename = file.path(output_plot_dir, "violin_plot_round_one_acgh_hgsa.png"),
  plot = violin_plot,
  width = 10,
  height = 5,
  dpi = 300,
  bg = "white"
)






