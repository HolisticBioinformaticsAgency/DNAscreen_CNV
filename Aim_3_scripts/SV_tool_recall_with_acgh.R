library(tidyverse)
library(scales)
library(ggplot2)
library(patchwork)

# ── Read ──────────────────────────────────────────────────────────────────────
df <- read_tsv("master_cnv_comparison_v3.tsv", show_col_types = FALSE)

parse_tool <- function(x) case_when(
  grepl("Sniffles", x) ~ "Sniffles2",
  grepl("CuteSV",   x) ~ "CuteSV",
  grepl("SVIM",     x) ~ "SVIM",
  grepl("DeBreak",  x) ~ "DeBreak"
)
parse_platform <- function(x) if_else(grepl("ONT", x), "ONT", "PacBio")

# ── Concordance columns ───────────────────────────────────────────────────────
concord_cols <- names(df)[grepl("concordant with aCGH", names(df))]

# ── Reshape + filter ──────────────────────────────────────────────────────────
combined <- df %>%
  filter(!`Sample ID` %in% c("Copy Neutral Male")) %>%
  select(Barcode, `Sample ID`, `CNV Type`, all_of(concord_cols)) %>%
  pivot_longer(cols      = all_of(concord_cols),
               names_to  = "col_name",
               values_to = "concordance") %>%
  mutate(
    tool     = parse_tool(col_name),
    platform = parse_platform(col_name)
  ) %>%
  # Exclude: tool not run (N/A), no aCGH reference, inconclusive aCGH
  filter(
    !is.na(concordance),
    !concordance %in% c("N/A", "No aCGH data", "Inconclusive aCGH data")
  ) %>%
  mutate(recalled = concordance == "YES")

# ── Recall summary ────────────────────────────────────────────────────────────
recall_summary <- combined %>%
  group_by(tool, platform) %>%
  summarise(
    n_total    = n(),
    n_recalled = sum(recalled, na.rm = TRUE),
    recall     = n_recalled / n_total,
    .groups    = "drop"
  ) %>%
  mutate(
    platform = factor(platform, levels = c("ONT", "PacBio")),
    label    = paste0(n_recalled, "/", n_total)
  )

# ── Print tables ──────────────────────────────────────────────────────────────
cat("\n── Recall summary (tool × platform) ──────────────────────────────────\n")
recall_summary %>%
  mutate(recall_pct = percent(recall, accuracy = 1)) %>%
  arrange(platform, desc(recall)) %>%
  select(platform, tool, n_recalled, n_total, recall_pct) %>%
  print(n = Inf)

cat("\n── Per-barcode recall detail ──────────────────────────────────────────\n")
combined %>%
  group_by(Barcode, `Sample ID`, `CNV Type`, tool, platform) %>%
  summarise(recalled = any(recalled, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from  = c(tool, platform),
              values_from = recalled,
              names_glue  = "{tool} ({platform})") %>%
  print(n = Inf)

# ── Plot builder ──────────────────────────────────────────────────────────────
tool_colours <- c(
  "Sniffles2" = "#1f77b4",
  "CuteSV"   = "#e6550d",
  "SVIM"     = "#31a354",
  "DeBreak"  = "#756bb1"
)

make_recall_plot <- function(data, plat, show_y_axis = TRUE) {
  d <- data %>%
    filter(platform == plat) %>%
    mutate(tool = fct_reorder(tool, recall, .desc = TRUE))
  
  p <- ggplot(d, aes(x = tool, y = recall, fill = tool)) +
    
    geom_hline(yintercept = 1, linetype = "dotted",
               colour = "grey50", linewidth = 0.5) +
    
    geom_col(width = 0.62, colour = "grey30", linewidth = 0.3) +
    
    geom_text(aes(label = label),
              vjust = -0.5, size = 3.5, fontface = "bold") +
    
    scale_fill_manual(values = tool_colours, guide = "none") +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1.12),
      breaks = seq(0, 1, 0.25)
    ) +
    
    labs(title = plat, x = NULL,
         y = if (show_y_axis) "Recall" else NULL) +
    
    theme_bw(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text.x        = element_text(size = 10),
      axis.text.y        = element_text(size = 9),
      axis.title.y       = element_text(size = 10),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  if (!show_y_axis) {
    p <- p + theme(axis.text.y  = element_blank(),
                   axis.ticks.y = element_blank())
  }
  p
}

p_ont    <- make_recall_plot(recall_summary, "ONT",    show_y_axis = TRUE)
p_pacbio <- make_recall_plot(recall_summary, "PacBio", show_y_axis = FALSE)

p_combined <- p_ont + p_pacbio +
  plot_annotation(
    title    = "CNV recall by SV caller and sequencing platform",
    subtitle = "Recall = concordant with aCGH.  Labels = recalled / total callable samples.",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

print(p_combined)