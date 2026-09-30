#!/usr/bin/env Rscript
# ============================================================
# STAGE 1: IDENTIFY AND REMOVE OUTLIER BAMs
# Port of automated_outlier_detection.R:
#   RPKM -> sample-sample correlation -> MDS (k = 2) -> DBSCAN
#   (eps from the elbow of the sorted k-NN distance curve)
#   DBSCAN noise points (cluster 0) are outliers.
# Samples with zero reads / constant coverage cannot be correlated and
# are flagged as outliers directly.
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(ggrepel)
  library(dbscan)
  library(plotly)
  library(htmlwidgets)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts", help = "counts RData from ReadInBams.R"),
  make_option("--run",    help = "run ID"),
  make_option("--minpts", type = "integer", default = 5L),
  make_option("--utils",  help = "path to framework_utils.R")
)))
source(opt$utils)

run    <- opt$run
counts <- load_counts(opt$counts)
raw    <- counts$raw                       # exons x samples
samples <- colnames(raw)

# ── RPKM (per sample) ───────────────────────────────────────────────────────
feature_kb <- (counts$anno$end - counts$anno$start + 1) / 1000
totals     <- colSums(raw)
rpkm       <- sweep(raw, 2, pmax(totals, 1), "/") / feature_kb * 1e6

usable <- totals > 0 & apply(rpkm, 2, function(x) sd(x) > 0)

status <- tibble(
  run         = run,
  sample      = samples,
  bam         = unname(counts$bams[samples]),
  total_reads = unname(totals),
  max_correlation = NA_real_,
  mds_x       = NA_real_,
  mds_y       = NA_real_,
  cluster     = NA_integer_,
  is_outlier  = !usable,
  outlier_reason = if_else(usable, NA_character_, "zero_or_constant_coverage")
)

if (sum(usable) >= opt$minpts + 1) {
  cor_matrix <- cor(rpkm[, usable, drop = FALSE])
  cor_no_diag <- cor_matrix
  diag(cor_no_diag) <- NA
  saveRDS(cor_no_diag, sprintf("%s_cor_matrix_no_diag.rds", run))

  max_cor <- apply(cor_no_diag, 2, max, na.rm = TRUE)

  mds <- cmdscale(1 - cor_matrix, k = 2)

  minPts <- opt$minpts
  k_dist <- apply(as.matrix(dist(mds)), 1, function(row) sort(row)[minPts])
  k_sorted <- sort(k_dist)
  elbow_index <- which.max(diff(diff(k_sorted))) + 1
  epsilon <- k_sorted[elbow_index]
  clustering <- dbscan(mds, eps = epsilon, minPts = minPts)
  message(sprintf("%s: DBSCAN eps = %.6f, minPts = %d, noise points = %d",
                  run, epsilon, minPts, sum(clustering$cluster == 0)))

  idx <- match(colnames(cor_matrix), status$sample)
  status$max_correlation[idx] <- max_cor
  status$mds_x[idx]   <- mds[, 1]
  status$mds_y[idx]   <- mds[, 2]
  status$cluster[idx] <- clustering$cluster
  status$is_outlier[idx]     <- clustering$cluster == 0
  status$outlier_reason[idx] <- if_else(clustering$cluster == 0, "dbscan_noise", NA_character_)

  # ── Histogram of max correlation values ──────────────────────────────────
  hist_p <- ggplot(tibble(max_correlation = max_cor), aes(x = max_correlation)) +
    geom_histogram(bins = 10, fill = "blue", color = "black") +
    labs(title = paste0("Histogram of Correlation Values for ", run),
         x = "Correlation", y = "Frequency") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white"))
  ggsave(sprintf("%s_correlation_histogram.png", run), hist_p, width = 8, height = 6)

  # ── MDS plot ──────────────────────────────────────────────────────────────
  plot_data <- status %>%
    filter(!is.na(cluster)) %>%
    mutate(Sample_ID = str_extract(sample, "DNS\\d+") %>% coalesce(sample),
           Cluster = factor(cluster))
  outliers <- plot_data %>% filter(cluster == 0)

  p <- ggplot(plot_data, aes(x = mds_x, y = mds_y, label = Sample_ID, color = Cluster)) +
    geom_point(size = 3) +
    geom_text_repel(data = outliers, aes(label = Sample_ID), size = 3,
                    segment.size = 0.2, segment.color = "black") +
    labs(title = paste0("Sample-sample correlation distance for DNA Screen ", run),
         x = "Dimension 1", y = "Dimension 2") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.background = element_rect(fill = "white", colour = NA)) +
    scale_color_discrete(name = "Cluster")
  ggsave(sprintf("%s_MDS_plot.png", run), p, width = 9, height = 7, dpi = 150)
  save_widget_html(ggplotly(p), sprintf("%s_MDS_plot.html", run))
} else {
  warning(sprintf("%s: only %d usable samples (< minPts + 1); DBSCAN skipped, no DBSCAN outliers called",
                  run, sum(usable)))
}

write_tsv_safe(status, sprintf("%s_outliers.tsv", run))
writeLines(status$bam[!status$is_outlier], sprintf("%s_retained_bams.txt", run))
write_tsv_safe(status %>% dplyr::select(run, sample, mds_x, mds_y, cluster, is_outlier),
               sprintf("%s_mds_coords.tsv", run))

message(sprintf("%s: %d BAMs, %d outliers, %d retained",
                run, nrow(status), sum(status$is_outlier), sum(!status$is_outlier)))
