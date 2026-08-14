library(tidyverse)
library(scales)
library(ggplot2)
library(ggtext)

MISS_THRESHOLD <- 2000
X_CAP          <- 2500

# ── Read ──────────────────────────────────────────────────────────────────────
df <- read_tsv("master_cnv_comparison_v3.tsv", show_col_types = FALSE) %>%
  mutate(cnv_size       = `True End` - `True Start`,
         cnv_size_label = paste0(format(cnv_size, big.mark = ","), " bp"))

# ── Desired barcode order (replicates grouped together) ───────────────────────
barcode_order <- c(
  "bc1002",                   # DNS006443
  "bc1008", "bc1101",         # DNS006916
  "bc1102",                   # DNS006946
  "bc1104",                   # DNS008387
  "bc1005",                   # DNS000601
  "bc1004", "bc1098",         # DNS004024
  "bc1003", "bc1100",         # DNS004937
  "bc1097",                   # DNS006683
  "bc1007"                    # DNS010349
)

# ── Helper: assign Rep. label based on bc number ──────────────────────────────
get_rep_label <- function(barcode) {
  bc_num <- as.integer(gsub("bc", "", barcode))
  case_when(
    bc_num >= 1001 & bc_num <= 1008 ~ "Rep. 1",
    bc_num >= 1097 & bc_num <= 1104 ~ "Rep. 2",
    TRUE ~ NA_character_
  )
}

# ── Build barcode_label: show Rep. only when Sample ID appears more than once ─
build_barcode_label <- function(d) {
  rep_lookup <- d %>%
    distinct(Barcode, `Sample ID`) %>%
    group_by(`Sample ID`) %>%
    mutate(n_bc = n()) %>%
    ungroup() %>%
    mutate(rep_label = if_else(n_bc > 1, get_rep_label(Barcode), NA_character_))
  
  d %>%
    left_join(rep_lookup, by = c("Barcode", "Sample ID")) %>%
    mutate(
      sample_id_line = if_else(
        !is.na(rep_label),
        paste0(`Sample ID`, " (", rep_label, ")"),
        `Sample ID`
      ),
      barcode_label = paste0(sample_id_line, "\n", Gene, " ", `CNV Type`, "\n(", cnv_size_label, ")")
    ) %>%
    select(-n_bc, -rep_label, -sample_id_line)
}

# ── "Did not call" from Coordinates columns ───────────────────────────────────
coord_cols <- names(df)[grepl("Coordinates on", names(df))]

did_not_call_raw <- df %>%
  filter(!`Sample ID` %in% c("Copy Neutral Male"),
         !Barcode %in% c("bc1006", "bc1103")) %>%
  select(Barcode, `Sample ID`, Gene, `CNV Type`, cnv_size_label, all_of(coord_cols)) %>%
  pivot_longer(cols = all_of(coord_cols),
               names_to = "col_name", values_to = "coord_val") %>%
  mutate(
    tool     = case_when(
      grepl("Sniffles", col_name) ~ "Sniffles",
      grepl("CuteSV",   col_name) ~ "CuteSV",
      grepl("SVIM",     col_name) ~ "SVIM"
    ),
    platform = if_else(grepl(" ONT$", col_name), "ONT", "PacBio")
  ) %>%
  filter(!is.na(tool),
         grepl("^did not call", trimws(coord_val), ignore.case = TRUE))

# ── Offsets ───────────────────────────────────────────────────────────────────
offset_cols <- names(df)[grepl("Offset of (Start|End)", names(df))]

long <- df %>%
  filter(!`Sample ID` %in% c("Copy Neutral Male"),
         !Barcode %in% c("bc1006", "bc1103")) %>%
  select(Barcode, `Sample ID`, Gene, `CNV Type`, cnv_size_label, all_of(offset_cols)) %>%
  pivot_longer(cols = all_of(offset_cols),
               names_to = "col_name", values_to = "offset_raw") %>%
  mutate(
    tool       = case_when(
      grepl("Sniffles", col_name) ~ "Sniffles",
      grepl("CuteSV",   col_name) ~ "CuteSV",
      grepl("SVIM",     col_name) ~ "SVIM"
    ),
    platform   = if_else(grepl("ONT", col_name), "ONT", "PacBio"),
    coord_type = if_else(grepl("Start coord", col_name), "start", "end"),
    offset     = suppressWarnings(as.numeric(offset_raw))
  ) %>%
  filter(!is.na(offset), !is.na(tool))

# ── Factor levels ─────────────────────────────────────────────────────────────
lvl_tool     <- c("Sniffles", "CuteSV", "SVIM")
lvl_platform <- c("ONT", "PacBio")
lvl_toolplat <- c(
  "Sniffles\n(ONT)", "Sniffles\n(PacBio)",
  "CuteSV\n(ONT)",   "CuteSV\n(PacBio)",
  "SVIM\n(ONT)",     "SVIM\n(PacBio)"
)

apply_factors <- function(d) {
  d %>%
    build_barcode_label() %>%
    mutate(
      tool      = factor(tool,     levels = lvl_tool),
      platform  = factor(platform, levels = lvl_platform),
      tool_plat = factor(paste0(tool, "\n(", platform, ")"), levels = lvl_toolplat)
    )
}

# ── Helper: build ordered barcode_label levels from barcode_order ─────────────
make_label_levels <- function(d, barcode_order) {
  ref <- d %>%
    distinct(Barcode, barcode_label) %>%
    mutate(Barcode = factor(Barcode, levels = barcode_order)) %>%
    arrange(Barcode)
  ref$barcode_label
}

# ── Combined offset ───────────────────────────────────────────────────────────
combined <- long %>%
  pivot_wider(
    id_cols    = c(Barcode, `Sample ID`, Gene, `CNV Type`, cnv_size_label, tool, platform),
    names_from = coord_type, values_from = offset
  ) %>%
  mutate(
    total_abs_offset = abs(start) + abs(end),
    is_miss          = total_abs_offset > MISS_THRESHOLD,
    offset_plot      = pmin(total_abs_offset, X_CAP)
  ) %>%
  apply_factors()

dnc_df <- did_not_call_raw %>%
  select(Barcode, `Sample ID`, Gene, `CNV Type`, cnv_size_label, tool, platform) %>%
  distinct() %>%
  apply_factors()

# ── Derive ordered barcode_label levels ───────────────────────────────────────
label_levels <- make_label_levels(combined, barcode_order)

# Re-apply ordered factor to barcode_label in all downstream data
reorder_label <- function(d) {
  d %>% mutate(barcode_label = factor(barcode_label, levels = label_levels))
}

combined    <- reorder_label(combined)
dnc_df      <- reorder_label(dnc_df)

# ── Valid calls ───────────────────────────────────────────────────────────────
valid_calls <- combined %>%
  filter(!is_miss) %>%
  mutate(plot_shape = as.character(platform), x_pos = offset_plot)

# ── No-valid-call rows (misses + did-not-call) → X marks at x = 0 ────────────
no_valid_calls <- bind_rows(
  combined %>%
    filter(is_miss) %>%
    mutate(plot_shape = "No valid call", x_pos = 0),
  dnc_df %>%
    mutate(offset_plot = 0, is_miss = FALSE,
           plot_shape = "No valid call", x_pos = 0)
) %>%
  distinct(Barcode, `Sample ID`, tool_plat, barcode_label,
           plot_shape, x_pos, .keep_all = TRUE)

all_points <- bind_rows(valid_calls, no_valid_calls) %>%
  mutate(
    plot_shape    = factor(plot_shape, levels = c("ONT", "PacBio", "No valid call")),
    barcode_label = factor(barcode_label, levels = label_levels)
  )

# ── Soft-clip dotted line data ────────────────────────────────────────────────
sc_cols_ont    <- c("Soft-clip Left clip pos offset vs Real BP Start (ONT)",
                    "Soft-clip Right clip pos offset vs Real BP End (ONT)")
sc_cols_pacbio <- c("Soft-clip Left clip pos offset vs Real BP Start (PACBIO)",
                    "Soft-clip Right clip pos offset vs Real BP End (PACBIO)")

softclip_all <- df %>%
  filter(!`Sample ID` %in% c("Copy Neutral Male"),
         !Barcode %in% c("bc1006", "bc1103")) %>%
  select(Barcode, `Sample ID`, Gene, `CNV Type`, cnv_size_label,
         sc_ont_left  = all_of(sc_cols_ont[1]),
         sc_ont_right = all_of(sc_cols_ont[2]),
         sc_pb_left   = all_of(sc_cols_pacbio[1]),
         sc_pb_right  = all_of(sc_cols_pacbio[2])) %>%
  pivot_longer(
    cols          = starts_with("sc_"),
    names_to      = c("platform_raw", "side"),
    names_pattern = "sc_(ont|pb)_(left|right)"
  ) %>%
  mutate(
    platform  = if_else(platform_raw == "ont", "ONT", "PacBio"),
    value_num = suppressWarnings(as.numeric(value))
  ) %>%
  group_by(Barcode, `Sample ID`, Gene, `CNV Type`, cnv_size_label, platform) %>%
  summarise(
    total_sc_abs = sum(abs(value_num), na.rm = TRUE),
    n_valid      = sum(!is.na(value_num)),
    .groups = "drop"
  ) %>%
  filter(n_valid == 2) %>%
  build_barcode_label() %>%                          # ← use helper here too
  mutate(
    barcode_label = factor(barcode_label, levels = label_levels),
    sc_x = pmin(total_sc_abs, X_CAP)
  )

# ── Scales ────────────────────────────────────────────────────────────────────
tool_colours <- c("Sniffles" = "#1f77b4", "CuteSV" = "#e6550d", "SVIM" = "#31a354")

shape_values <- c("ONT" = 16, "PacBio" = 17, "No valid call" = 4)
shape_labels <- c(
  "ONT"           = "ONT",
  "PacBio"        = "PacBio",
  "No valid call" = "✕  No valid call\n(not found in region, or offset > 2,000 bp)"
)

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot() +
  
  geom_vline(
    data      = softclip_all %>% filter(platform == "ONT"),
    aes(xintercept = sc_x),
    colour    = "#add8e6",
    linetype  = "dotted",
    linewidth = 0.9
  ) +
  geom_vline(
    data      = softclip_all %>% filter(platform == "PacBio"),
    aes(xintercept = sc_x),
    colour    = "#4b0082",
    linetype  = "dotted",
    linewidth = 0.9
  ) +
  
  geom_segment(
    data = valid_calls,
    aes(x = 0, xend = x_pos, y = tool_plat, yend = tool_plat, colour = tool),
    linewidth = 0.55, alpha = 0.55, show.legend = FALSE
  ) +
  
  geom_point(
    data = all_points %>% filter(plot_shape != "No valid call"),
    aes(x = x_pos, y = tool_plat, colour = tool, shape = plot_shape),
    size = 3, stroke = 0.6
  ) +
  
  geom_point(
    data = all_points %>% filter(plot_shape == "No valid call"),
    aes(x = x_pos, y = tool_plat, colour = tool, shape = plot_shape),
    size = 3.5, stroke = 1.3
  ) +
  
  facet_wrap(~ barcode_label, nrow = 2) +
  
  scale_colour_manual(values = tool_colours, name = "Tool") +
  scale_shape_manual(
    values = shape_values,
    labels = shape_labels,
    name   = "Platform / Call status"
  ) +
  scale_x_continuous(
    limits = c(0, X_CAP),
    breaks = c(0, 500, 1000, 1500, 2000, 2500),
    labels = c("0", "500", "1,000", "1,500", "2,000", ">2,000\n(miss)"),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  
  labs(
    title    = "Total absolute offset of SV caller and consensus soft-clipped coordinates relative to manually-resolved breakpoints",
    subtitle = "Dotted lines represent offset of consensus soft-clipped coordinates (light blue = ONT,  dark purple = PacBio)",
    x        = "Total absolute offset (bp)",
    y        = NULL
  ) +
  
  theme_bw(base_size = 10) +
  theme(
    strip.text         = element_text(size = 7, face = "bold"),
    strip.background   = element_rect(fill = "grey93", colour = "grey70"),
    axis.text.y        = element_text(size = 7.5),
    axis.text.x        = element_text(size = 7, angle = 30, hjust = 1),
    axis.title.x       = element_text(size = 9),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold", size = 9),
    legend.text        = element_text(size = 8),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey88"),
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 8, colour = "grey40"),
    plot.margin        = margin(10, 14, 10, 14)
  )

print(p)