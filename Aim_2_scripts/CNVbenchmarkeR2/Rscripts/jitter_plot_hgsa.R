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
# zscore_combined <- zscore_combined %>%
#   filter(acgh_status != "tbc")

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
empty_tools <- tibble(
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

zscore_combined <- bind_rows(zscore_combined, empty_tools)

# 3️⃣ Apply the consistent factor ordering
zscore_combined$Tool <- factor(zscore_combined$Tool, levels = tool_order)


# Plot all tools in one figure
combined_plot <- ggplot(zscore_combined, aes(x = Tool, y = zscore)) +
  # First: plot all points except the two special categories
  geom_jitter(
    data = subset(zscore_combined, 
                  !Highlight %in% c("BRCA1 deletion concordant with aCGH",
                                    # "PCSK9 duplication concordant with aCGH",
                                    "APOB duplication discordant with aCGH")),
                                    # "PCSK9 duplication discordant with aCGH"
    aes(color = Highlight, shape = Highlight),
    width = 0.2, alpha = 0.6, size = 1.5
  ) +
  # Second: plot the two special categories on top
  geom_jitter(
    data = subset(zscore_combined, 
                  Highlight %in% c("BRCA1 deletion concordant with aCGH",
                                   # "PCSK9 duplication concordant with aCGH",
                                   "APOB duplication discordant with aCGH")),
                                   # "PCSK9 duplication discordant with aCGH"
    aes(color = Highlight, shape = Highlight),
    width = 0.2, size = 2.5 # slightly bigger so they stand out
  ) +
  scale_color_manual(values = c(
    "Concordant with aCGH" = "green",
    "Discordant with aCGH" = "red",
    "BRCA1 deletion concordant with aCGH" = "#7570B3",
    # "PCSK9 duplication concordant with aCGH" = "#FF7F00",
    "APOB duplication discordant with aCGH" = "#E6AB02",
    # "PCSK9 duplication discordant with aCGH" = "#377EB8",
    "From a Sample with tbc aCGH status" = "#666666"
  )) +
  scale_shape_manual(values = c(
    "Concordant with aCGH" = 17,
    "Discordant with aCGH" = 17,
    "BRCA1 deletion concordant with aCGH" = 16,
    "PCSK9 duplication concordant with aCGH" = 16,
    "APOB duplication discordant with aCGH" = 16,
    "PCSK9 duplication discordant with aCGH" = 16,
    "From a Sample with tbc aCGH status" = 16
  )) +
  theme_minimal(base_size = 12) +
  labs(
    title = NULL,
    x = "Tool",
    y = "Z-scaled Quality Score",
    color = "CNV",
    shape = "CNV"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",        # move to bottom
    legend.direction = "horizontal",   # horizontal order
    legend.box = "horizontal"          # align horizontally
  )

combined_plot

# Save to file
# ggsave(
#   filename = file.path(output_plot_dir, "zscore_combined_tools_round_one_acgh_hgsa.png"),
#   plot = combined_plot,
#   width = 10,
#   height = 5,
#   dpi = 300,
#   bg = "white"
# )


###### PLOTLY ######
combined_plot <- ggplot(zscore_combined,
                        aes(x = Tool, y = zscore, color = Highlight,
                            text = paste("Sample:", Sample,
                                         "<br>Gene:", Gene,
                                         "<br>Coordinates:", paste0("chr", Chromosome, ":", Start, "-", End),
                                         "<br>Raw score:", round(raw_score, 3),
                                         "<br>Z-score:", round(zscore, 2),
                                         "<br>Status:", Highlight))) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Concordant with aCGH" = "green",
    "Discordant with aCGH" = "red",
    "BRCA1 deletion concordant with aCGH" = "#7570B3",
    "PCSK9 duplication concordant with aCGH" = "#FF7F00",
    "APOB duplication discordant with aCGH" = "#E6AB02",
    "PCSK9 duplication discordant with aCGH" = "#377EB8",
    "From a Sample with tbc aCGH status" = "#666666")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Z-scaled Quality Scores per Tool",
    x = "Tool",
    y = "Z-scaled Quality Score",
    color = "CNV"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

# Use 'tooltip = "text"' to activate the custom hover text
interactive_plot <- ggplotly(combined_plot, tooltip = "text")
interactive_plot

saveWidget(interactive_plot, file = file.path(output_plot_dir, "zscore_combined_tools_round_one_acgh_with_unnconfirmed.html"), selfcontained = TRUE)


