library(dplyr)
library(ggplot2)
library(tidyr)
library(IRanges)  # for overlap detection
library(stringr)
library(GenomicRanges)
library(patchwork)


# Step 1: initiate the matrix NOTE THAT THIS IS BASED ON VARSEQ CNVS FIRST, does not include other CNVs found uniquely
# by aCGH or other tools... (there aren't any nonvs-cnv calls from other tools that are validated by aCGH either...)

# Step 1: Initialize result matrix
# Start with CNVs sent for aCGH ONLY! 
# Validated_res are CNVs that are the 92 candidates selected for aCGH!
# cnv_from_acgh_sample <- vscnv_reorder_perrun %>%
#   filter(Sample %in% unique(validated_res$sample)) %>%
#   select(Sample, Gene, CNV_Type, Chromosome, Start, End, coordinates)

cnv_from_acgh_sample <- validated_res %>%
  mutate(CNV_Type = case_when(
    (Type.of.CNV == "Het Deletion") ~ "deletion",
    (Type.of.CNV == "Duplicate") ~ "duplication",
    TRUE ~ NA_character_
  )) %>%
  rename(Sample = sample, 
         Gene = Gene.Names,
         Chromosome = VS_Call_chr,
         Start = VS_Call_start,
         End = VS_Call_end,
         coordinates = coordinates) %>%
  select(Sample, Gene, CNV_Type, Chromosome, Start, End, coordinates)


# Step 2: Initialize result matrix with tool columns and acgh_confirmed
cnv_result_matrix <- cnv_from_acgh_sample

# Add one column per tool, initialized to NA
for (tool_name in names(harmonised_list_with_acgh)) {
  cnv_result_matrix[[tool_name]] <- NA
}

# Add acgh_confirmed column
cnv_result_matrix$acgh_confirmed <- NA

# Step 3: Loop through each CNV
for (i in seq_len(nrow(cnv_result_matrix))) {
  cnv_row <- cnv_result_matrix[i, ]
  ref_range <- IRanges(start = cnv_row$Start, end = cnv_row$End)
  confirmed_by_any_tool <- FALSE
  
  for (tool_name in names(harmonised_list_with_acgh)) {
    tool_data <- harmonised_list_with_acgh[[tool_name]]$data
    
    # Filter tool data by Sample and Chromosome
    tool_calls <- tool_data %>%
      filter(Sample == cnv_row$Sample,
             Chromosome == cnv_row$Chromosome)
    
    if (nrow(tool_calls) > 0) {
      tool_ranges <- IRanges(start = tool_calls$Start, end = tool_calls$End)
      hits <- findOverlaps(ref_range, tool_ranges)
      
      if (length(hits) > 0) {
        # Get overlapping tool call(s)
        overlap_idx <- subjectHits(hits)
        acgh_statuses <- tool_calls$acgh_status[overlap_idx]
        
        # Check if any are confirmed
        if ("concordant" %in% acgh_statuses) {
          confirmed_by_any_tool <- TRUE
          cnv_result_matrix[i, tool_name] <- "TP_by_tool"
        } else if ("discordant" %in% acgh_statuses) {
          cnv_result_matrix[i, tool_name] <- "FP_by_tool"
        }
        else if ("tbc" %in% acgh_statuses) {
          cnv_result_matrix[i, tool_name] <- "tbc" # caller calls it but sample acgh status tbc...
          confirmed_by_any_tool <- "tbc" # as long as one of the cnv caller has tbc, which means this sample status is not confirmed yet...
        }
        else if ("inconclusive" %in% acgh_statuses) {
          cnv_result_matrix[i, tool_name] <- "inconclusive" # caller calls it but sample acgh status tbc...
          confirmed_by_any_tool <- "inconclusive" # as long as one of the cnv caller has tbc, which means this sample status is not confirmed yet...
        }
      } else {
        cnv_result_matrix[i, tool_name] <- "Not_Called"
      }
    } else {
      cnv_result_matrix[i, tool_name] <- "Not_Called"
    }
  }
  
  # Set acgh_confirmed for this CNV row
  cnv_result_matrix$acgh_confirmed[i] <- confirmed_by_any_tool
}


# Clean up tool naming
cnv_result_matrix <- cnv_result_matrix %>%
  select(-"VS-CNV_val")

colnames(cnv_result_matrix)[colnames(cnv_result_matrix) == "VS-CNV_perrun"] <- "VS-CNV"

##### FILTER OUR TBC or inconclusive CNVs for Aims 1 and 2...
cnv_result_matrix <- cnv_result_matrix %>%
  filter(!acgh_confirmed %in% c("tbc", "inconclusive"))

# Step 3: Simplify Sample names
cnv_result_matrix <- cnv_result_matrix %>%
  mutate(Sample = sub(".*(DNS\\d+).*", "\\1", Sample))

# Step 4: Format aCGH column
cnv_result_matrix <- cnv_result_matrix %>%
  mutate(aCGH = as.character(acgh_confirmed)) %>%
  select(-acgh_confirmed)

# Step 5: Get tool columns (all columns beyond fixed ones)
tool_cols <- colnames(cnv_result_matrix)[!(colnames(cnv_result_matrix) %in%
                                             c("Sample", "Gene", "CNV_Type", "Chromosome", "Start", "End", "coordinates"))]

# Count tool TP calls
cnv_result_matrix <- cnv_result_matrix %>%
  mutate(ToolCalls = rowSums(across(all_of(tool_cols)) == "TP_by_tool", na.rm = TRUE))

# Step 6: Long format
cnv_call_long <- cnv_result_matrix %>%
  unite("CNV_ID", Sample, Gene, CNV_Type, Start, End, sep = " | ", remove = FALSE) %>%
  pivot_longer(
    cols = all_of(tool_cols),
    names_to = "Tool",
    values_to = "CallStatus"
  )

# Step 7: Assign fill status
cnv_call_long <- cnv_call_long %>%
  mutate(
    FillStatus = case_when(
      CallStatus == "TP_by_tool" ~ "TP_by_tool",
      CallStatus == "FP_by_tool" ~ "FP_by_tool",
      Tool == "aCGH" & CallStatus == "TRUE" ~ "aCGH_TP",
      Tool == "aCGH" & CallStatus == "FALSE" ~ "aCGH_FP",
      # CallStatus == "tbc" ~ "tbc",
      CallStatus == "Not_Called" ~ "Not_called"
    )
  ) %>%
  mutate(FillStatus = if_else(Tool == "VS-CNV", "Candidates", FillStatus))
  

# Step 8: Order CNV_IDs by TP support
cnv_id_order <- cnv_result_matrix %>%
  unite("CNV_ID", Sample, Gene, CNV_Type, Start, End, sep = " | ", remove = FALSE) %>%
  # arrange(desc(ToolCalls)) %>%
  pull(CNV_ID)

cnv_call_long$CNV_ID <- factor(cnv_call_long$CNV_ID, levels = cnv_id_order)

# Step 9: Order tools (put aCGH last)
# Count TP_by_tool calls per tool
tool_tp_counts <- cnv_call_long %>%
  filter(CallStatus == "TP_by_tool") %>%
  count(Tool, name = "TP_count") %>%
  arrange(desc(TP_count))

# Arrange tools: VS-CNV first, aCGH second, rest by TP_count
other_tools <- tool_tp_counts %>%
  filter(!Tool %in% c("VS-CNV", "aCGH")) %>%
  arrange(desc(TP_count)) %>%
  pull(Tool)

tool_order <- c("VS-CNV", "Spacer1", "aCGH", "Spacer2", other_tools)

# Get unique CNV_IDs
unique_ids <- unique(cnv_call_long$CNV_ID)

# Add blank Spacer rows to cnv_call_long
spacer1_rows <- data.frame(
  CNV_ID = unique_ids,
  Tool = "Spacer1",
  CallStatus = NA,
  FillStatus = NA
)

spacer2_rows <- data.frame(
  CNV_ID = unique_ids,
  Tool = "Spacer2",
  CallStatus = NA,
  FillStatus = NA
)

cnv_call_long <- bind_rows(cnv_call_long, spacer1_rows, spacer2_rows)

# Apply ordering
cnv_call_long$Tool <- factor(cnv_call_long$Tool, levels = tool_order)


# X bounds: from first to last Tool
x_min <- 0.5
x_max <- length(levels(cnv_call_long$Tool)) + 0.5

# Add box to the plot
concordance_plot <- ggplot(cnv_call_long, aes(x = Tool, y = CNV_ID, fill = FillStatus)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c(
      "TP_by_tool" = "#1b9e77",
      "FP_by_tool" = "red",
      "aCGH_TP"    = "#1b9e77",
      "aCGH_FP"    = "red",
      "Not_called" = "grey",
      "Candidates" = "grey",
      "Spacer1"     = "white",
      "Spacer2"     = "white"
    ),
    na.translate = FALSE
  ) +
  scale_y_discrete(limits = rev(levels(cnv_call_long$CNV_ID))) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = NULL,
    legend.position = "none"
  ) +
  labs(x = NULL, y = NULL) + 
  scale_x_discrete(labels = function(x) {
    ifelse(x %in% c("Spacer1", "Spacer2"), "", x)
  })


# Add dummy rows for spacers (if not already in the data)
spacer_rows <- data.frame(
  algorithm = c("Spacer1", "Spacer2", "aCGH"),
  value = 0,
  segment = "FP"
)

gene_fp_df <- bind_rows(gene_fp_df, spacer_rows)

gene_fp_df$algorithm <- factor(gene_fp_df$algorithm, levels = tool_order)

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
        axis.title.x = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


p_gene_inverted <- p_gene +
  scale_y_reverse(expand = expansion(mult = c(0, 0.1))) +  # invert y-axis
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_blank(),  # optional: hide x labels
    axis.ticks.x = element_blank(),
  )

# Adjust spacing with theme + reduce vertical gap
concordance_plot <- concordance_plot + theme(plot.margin = margin(5, 5, 5, 5))
p_gene_inverted <- p_gene_inverted + theme(plot.margin = margin(0, 5, 5, 5))  # Top margin 0

concordance_plot <- concordance_plot / p_gene_inverted +
  plot_layout(heights = c(3, 1), guides = "collect")

print(concordance_plot)



####### Top 20 CNVs from concordance plot #######
# STEP 1: Get top 20 CNVs (by frequency or any other logic)

# Get top 20 CNVs
cnv_id_order <- unique(cnv_call_long$CNV_ID)[1:21]

top20_data <- cnv_call_long %>%
  filter(CNV_ID %in% cnv_id_order) %>%
  mutate(CNV_ID = factor(CNV_ID, levels = rev(cnv_id_order)))  # REV for top at top

# Labels for top 20 CNVs
cnv_labels <- cnv_call_long %>%
  filter(CNV_ID %in% cnv_id_order) %>%
  distinct(CNV_ID, Gene, CNV_Type, Start, End) %>%
  mutate(
    CNV_Label = paste(Gene, CNV_Type, paste0(round((End - Start + 1)/1000, 2), "kbp"), sep = " | ")
  )

# Mapping CNV_ID to label
label_mapping <- setNames(cnv_labels$CNV_Label, cnv_labels$CNV_ID)

# Plot
concordance_plot_top20 <- top20_data %>%
  ggplot(aes(x = Tool, y = CNV_ID, fill = FillStatus)) +
  geom_tile(color = "white") +
  geom_text(
    data = top20_data %>%
      filter(Tool == "VS-CNV") %>%
      mutate(label_num = row_number()),
    aes(label = label_num),
    color = "black",
    size = 6
  ) +
  scale_fill_manual(
    values = c(
      "TP_by_tool" = "#1b9e77",
      "FP_by_tool" = "red",
      "aCGH_TP"    = "#1b9e77",
      "aCGH_FP"    = "red",
      "Not_called" = "grey",
      "Candidates" = "grey",
      "Spacer1"     = "white",
      "Spacer2"     = "white"
    ),
    labels = c(
      "TP_by_tool" = "Concordant with aCGH",
      "FP_by_tool" = "Discordant with aCGH",
      "aCGH_TP" = "Observed by aCGH",
      "aCGH_FP" = "Not observed by aCGH",
      "Not_called" = "Not Called",
      "Candidates" = "Candidates"
    ),
    na.translate = FALSE,
    name = NULL
  ) +
  scale_y_discrete(
    labels = label_mapping[1:21],
    position = "right"  # shift axis to the right
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y.left = element_blank(),         # hide left y-axis labels
    axis.ticks.y.left = element_blank(),        # hide left ticks
    axis.text.y.right = element_text(size = 16), # show right labels
    axis.ticks.y.right = element_line(),        # show right ticks
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 16)
  ) +
  labs(x = NULL, y = NULL) + scale_x_discrete(labels = function(x) {
    ifelse(x %in% c("Spacer1", "Spacer2"), "", x)
  })


print(concordance_plot_top20)

# # Save to file
# ggsave(
#   filename = file.path(dir, "summary/cnv_tool_concordance_hgsa.png"),
#   plot = concordance_plot,
#   width = 8,
#   height = 18,
#   dpi = 300,
#   bg = "white"
# )
# 
# # Save to file
# ggsave(
#   filename = file.path(dir, "summary/cnv_tool_concordance_top20_hgsa.png"),
#   plot = concordance_plot_top20,
#   width = 16,
#   height = 10,
#   dpi = 300,
#   bg = "white"
# )

