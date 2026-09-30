#!/usr/bin/env Rscript
# ============================================================
# Interactive MDS of all runs (port of correlation_matrix.R)
# Each run's MDS is computed within that run (coordinates are per run);
# click a run to highlight it, outliers drawn as crosses.
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(plotly)
  library(crosstalk)
  library(htmlwidgets)
})
opt <- parse_args(OptionParser(option_list = list(make_option("--utils"))))
source(opt$utils)

files <- list.files(".", pattern = "_mds_coords\\.tsv$")
coords <- bind_rows(lapply(files, function(f)
  suppressMessages(read_tsv(f, col_types = cols(.default = col_character()))))) %>%
  type_convert(guess_integer = TRUE) %>%
  filter(!is.na(mds_x))

if (nrow(coords) == 0) {
  message("No MDS coordinates to plot")
  quit(save = "no")
}

runs <- unique(coords$run)
runs <- runs[order(suppressWarnings(as.numeric(str_extract(runs, "\\d+"))))]
coords <- coords %>%
  mutate(DNAscreen_Run = factor(run, levels = runs),
         Status = if_else(is_outlier %in% TRUE, "outlier", "retained"),
         Sample_ID = sample)

shared <- SharedData$new(coords, key = ~DNAscreen_Run)
p <- plot_ly(
  shared, x = ~mds_x, y = ~mds_y, type = "scatter", mode = "markers",
  color = ~DNAscreen_Run,
  colors = colorRampPalette(c("lightblue", "darkblue"))(length(runs)),
  symbol = ~Status, symbols = c("circle", "x"),
  text = ~paste("Sample ID:", Sample_ID, "<br>Run:", DNAscreen_Run, "<br>", Status),
  hoverinfo = "text"
) %>%
  layout(title = "MDS Plot of Samples (Filter by Run)",
         xaxis = list(title = "Dimension 1"), yaxis = list(title = "Dimension 2")) %>%
  highlight(on = "plotly_click", off = "plotly_doubleclick", dynamic = TRUE, color = "red")

save_widget_html(p, "MDS_all_runs.html")
