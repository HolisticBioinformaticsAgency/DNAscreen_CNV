library(dplyr)
library(ggplot2)
library(tidyr)
library(IRanges)  # for overlap detection
library(stringr)
library(GenomicRanges)


# ### Including CNVs found by aCGH but not found by VS-CNV per run


# Step 1: initiate the matrix NOTE THAT THIS IS BASED ON VARSEQ CNVS FIRST, does not include other CNVs found uniquely
# by aCGH or other tools...


# Step 1: Initialize result matrix
# Step 1: Start with CNVs sent for aCGH
cnv_selected_for_acgh <- vscnv_reorder_perrun %>%
  filter(cnv_selected_for_acgh == TRUE) %>%
  dplyr::select(Sample, Gene, CNV_Type, Chromosome, Start, End, coordinates)

# Step 2: Initialize result matrix with tool columns and acgh_confirmed
cnv_result_matrix <- cnv_selected_for_acgh

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
          confirmed_by_any_tool <- NA # as long as one of the cnv caller has tbc, which means this sample status is not confirmed yet...
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
  dplyr::select(-"VS-CNV_val")

colnames(cnv_result_matrix)[colnames(cnv_result_matrix) == "VS-CNV_perrun"] <- "VS-CNV"

##### FILTER OUR TBC CNVs FOR NOW!
### For now showing it internally (suffixed output files "with_unconfirmed")
# cnv_result_matrix <- cnv_result_matrix %>%
#   filter(!is.na(acgh_confirmed))

# Step 3: Simplify Sample names
cnv_result_matrix <- cnv_result_matrix %>%
  mutate(Sample = sub(".*(DNS\\d+).*", "\\1", Sample))

# Step 4: Format aCGH column
cnv_result_matrix <- cnv_result_matrix %>%
  mutate(aCGH = as.character(acgh_confirmed)) %>%
  dplyr::select(-acgh_confirmed)
cnv_result_matrix$aCGH[is.na(cnv_result_matrix$aCGH)] <- "tbc"

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

# im using a sample list instead of cnv list as a shortcut....
sample_with_unconfirmed_cnv = c("DNSxxxxxx",
                                "DNSxxxxxx",
                                "DNSxxxxxx",
                                "DNSxxxxxx",
                                "DNSxxxxxx",
                                "DNSxxxxxx",
                                "DNSxxxxxx",
                                "DNSxxxxxx")

# Step 7: Assign fill status
cnv_call_long <- cnv_call_long %>%
  mutate(
    FillStatus = case_when(
      CallStatus == "tbc" & Sample %in% sample_with_unconfirmed_cnv ~ "From a tbc CNV",
      CallStatus == "tbc" ~ "From a Sample not on array yet",
      CallStatus == "TP_by_tool" ~ "TP_by_tool",
      CallStatus == "FP_by_tool" ~ "FP_by_tool",
      Tool == "aCGH" & CallStatus == "TRUE" ~ "aCGH_TP",
      Tool == "aCGH" & CallStatus == "FALSE" ~ "aCGH_FP",
      CallStatus == "Not_Called" ~ "Not_called"
    )
  )

# Step 8: Order CNV_IDs by TP support
cnv_id_order <- cnv_result_matrix %>%
  unite("CNV_ID", Sample, Gene, CNV_Type, Start, End, sep = " | ", remove = FALSE) %>%
  arrange(desc(ToolCalls)) %>%
  pull(CNV_ID)

cnv_call_long$CNV_ID <- factor(cnv_call_long$CNV_ID, levels = cnv_id_order)

# Top 20 CNVs
top20_ids <- unique(cnv_call_long$CNV_ID)[1:20]  # adjust based on your logic

# Step 9: Order tools (put aCGH last)
# Count TP_by_tool calls per tool
tool_tp_counts <- cnv_call_long %>%
  filter(CallStatus == "TP_by_tool") %>%
  dplyr::count(Tool, name = "TP_count") %>%
  arrange(desc(TP_count))

# Arrange tools: VS-CNV first, aCGH second, rest by TP_count
other_tools <- tool_tp_counts %>%
  filter(!Tool %in% c("VS-CNV", "aCGH")) %>%
  arrange(desc(TP_count)) %>%
  pull(Tool)

tool_order <- c("VS-CNV", "aCGH", other_tools)

# Apply ordering
cnv_call_long$Tool <- factor(cnv_call_long$Tool, levels = tool_order)


# Recalculate y-axis position for top 20
all_ids <- levels(cnv_call_long$CNV_ID)
y_positions <- which(all_ids %in% top20_ids)
y_min_box <- min(y_positions) - 0.5
y_max_box <- max(y_positions) + 0.5

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
      "aCGH_TP"    = "darkgreen",
      "aCGH_FP"    = "darkred",
      "From a tbc CNV" = "pink",
      "Not_called" = "grey",
      "From a Sample not on array yet" = "darkgrey"
    ),
    na.translate = FALSE
  ) +
  scale_y_discrete(limits = rev(levels(cnv_call_long$CNV_ID))) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = NULL)

# Get x-position of the 'aCGH' column
acgh_pos <- which(levels(cnv_call_long$Tool) == "aCGH")

# Define y-axis bounds (entire height of CNV_IDs)
y_min <- 0.5
y_max <- length(unique(cnv_call_long$CNV_ID)) + 0.5

# Add rectangle layer for the border
concordance_plot <- concordance_plot + annotate("rect",
             xmin = acgh_pos - 0.5, xmax = acgh_pos + 0.5,
             ymin = y_min, ymax = y_max,
             fill = NA, color = "black", size = 1.5)


# Display plot
print(concordance_plot)


####### Top 20 CNVs from concordance plot #######
# STEP 1: Get top 20 CNVs (by frequency or any other logic)
top20_cnvs <- cnv_call_long %>%
  dplyr::count(CNV_ID, sort = TRUE) %>%
  slice_head(n = 20) %>%
  pull(CNV_ID)


# Only unique CNVs (ignore sample)
cnv_labels <- cnv_call_long %>%
  filter(CNV_ID %in% top20_cnvs) %>%
  distinct(CNV_ID, Sample, Gene, CNV_Type, Start, End) %>%
  mutate(
    CNV_Label = paste(Sample, Gene, CNV_Type, paste0(End - Start + 1, "bp"), sep = " | ")
  )

cnv_id_order <- cnv_call_long %>%
  dplyr::count(CNV_ID, sort = TRUE) %>%
  slice_head(n = 20) %>%
  pull(CNV_ID)

# Filter for top 20 and set CNV_ID as factor with this order
top20_data <- cnv_call_long %>%
  filter(CNV_ID %in% cnv_id_order) %>%
  mutate(CNV_ID = factor(CNV_ID, levels = rev(cnv_id_order)))  # REV for top at top


# Mapping CNV_ID to labels
label_mapping <- setNames(cnv_labels$CNV_Label, cnv_labels$CNV_ID)

# Plot
concordance_plot_top20 <- top20_data %>%
  ggplot(aes(x = Tool, y = CNV_ID, fill = FillStatus)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c(
      "TP_by_tool" = "#1b9e77",
      "FP_by_tool" = "red",
      "aCGH_TP"    = "darkgreen",
      "aCGH_FP"    = "darkred",
      "From a tbc CNV" = "pink",
      "Not_called" = "grey",
      "From a Sample not on array yet" = "darkgrey"
    ),
    labels = c(
      "TP_by_tool" = "Concordant with aCGH",
      "FP_by_tool" = "Discordant with aCGH",
      "aCGH_TP" = "Observed by aCGH",
      "aCGH_FP" = "Not observed by aCGH",
      "Not_called" = "Not Called"
    ),
    na.translate = FALSE,
    name = NULL
  ) +
  scale_y_discrete(
    labels = label_mapping,
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
  labs(x = NULL, y = NULL)


# Get x-position of the 'aCGH' column
acgh_pos <- which(levels(cnv_call_long$Tool) == "aCGH")

# Define y-axis bounds (entire height of CNV_IDs)
y_min <- 0.5
y_max <- 20 + 0.5

# Add rectangle layer for the border
concordance_plot_top20 <- concordance_plot_top20 + annotate("rect",
                                                xmin = acgh_pos - 0.5, xmax = acgh_pos + 0.5,
                                                ymin = y_min, ymax = y_max,
                                                fill = NA, color = "black", size = 1.5)

print(concordance_plot_top20)

# # Save to file
# ggsave(
#   filename = file.path(dir, "summary/cnv_tool_concordance_with_unconfirmed.png"),
#   plot = concordance_plot,
#   width = 5,
#   height = 7,
#   dpi = 300,
#   bg = "white"
# )
# 
# # Save to file
# ggsave(
#   filename = file.path(dir, "summary/cnv_tool_concordance_top20_with_unconfirmed.png"),
#   plot = concordance_plot_top20,
#   width = 18,
#   height = 10,
#   dpi = 300,
#   bg = "white"
# )

# write.csv(cnv_result_matrix, file.path(dir, "summary/cnv_result_matrix.csv"))

