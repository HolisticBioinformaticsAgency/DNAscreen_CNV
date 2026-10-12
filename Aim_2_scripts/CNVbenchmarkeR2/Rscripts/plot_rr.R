library(ExomeDepth)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(tidyr)
library(readr)
library(GenomicRanges)
library(IRanges)

# ============================================================
# OUTPUT DIRECTORY
# ============================================================

output_dir <- "/path/to/project/Aim_2/Optimized_CNV_Detection_Pipeline/test"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# LOAD CACHED EXON DATA (gene track)
# ============================================================

exons_combined <- readRDS("/path/to/project/Aim_2/CNV-plot/cached_exon_data.rds")

clinically_relevant_transcripts <- list(
  BRCA1 = "NM_007294.4",
  BRCA2 = "NM_000059.4",
  PALB2 = "NM_024675.4",
  MSH2  = "NM_000251.3",
  MSH6  = "NM_000179.3",
  MLH1  = "NM_000249.4",
  APOB  = "NM_000384.3",
  PCSK9 = "NM_174936.4",
  LDLR  = "NM_000527.5"
)

# ============================================================
# MAIN PLOTTING FUNCTION
# ============================================================

plot_rr <- function(target_sample, gene_symbol, cnv_chr,
                          cnv_start, cnv_end, cnv_type,
                          dnascreen_run, output_dir,
                          exons_combined,
                          clinically_relevant_transcripts) {
  
  read_count_dir <- "/path/to/project/Aim_1/CNV_calling/exome_depth/RData"
  read_count_rdata <- sprintf("%s/run%d_bams_9genes_25bp.fix.sorted.RData",
                              read_count_dir, dnascreen_run)
  
  if (!file.exists(read_count_rdata)) {
    warning(sprintf("RData not found for run %d, skipping.", dnascreen_run))
    return(invisible(NULL))
  }
  
  load(read_count_rdata)  # loads my.count
  
  my.count.dafr <- as(my.count, "data.frame")
  
  sample_ids <- colnames(my.count.dafr)[6:ncol(my.count.dafr)]
  count_matrix <- my.count.dafr[, sample_ids] %>%
    mutate_all(as.numeric)
  
  # CPT normalisation only
  sample_totals <- colSums(count_matrix, na.rm = TRUE)
  normalized_counts <- sweep(count_matrix, 2, sample_totals, FUN = "/") * 1e3
  normalized_df <- as.data.frame(normalized_counts)
  colnames(normalized_df) <- sample_ids
  my.count.dafr[, sample_ids] <- normalized_df
  
  my.count.dafr$chromosome <- gsub("^chr", "", my.count.dafr$chromosome)
  my.count.dafr$chromosome <- as.numeric(my.count.dafr$chromosome)
  
  # decon_with_acgh Sample usually lacks suffix, so add it back
  target_sample_with_suffix <- paste0(target_sample, ".hq.sorted.marked.bam")
  target_sample_dots <- gsub("-", ".", target_sample_with_suffix)
  
  if (!target_sample_dots %in% sample_ids) {
    warning(sprintf("Sample %s not found in run %d, skipping.",
                    target_sample_dots, dnascreen_run))
    return(invisible(NULL))
  }
  
  gene_exon_coords <- exons_combined %>%
    filter(gene == gene_symbol,
           chromosome_name == as.character(cnv_chr))
  
  if (nrow(gene_exon_coords) == 0) {
    user_start <- cnv_start - 5000
    user_end   <- cnv_end + 5000
  } else {
    user_start <- min(gene_exon_coords$exon_chrom_start)
    user_end   <- max(gene_exon_coords$exon_chrom_end)
  }
  
  user_chr <- as.numeric(cnv_chr)
  
  filtered_counts <- my.count.dafr %>%
    dplyr::filter(chromosome == user_chr,
                  start >= user_start,
                  end <= user_end)
  
  if (nrow(filtered_counts) == 0) {
    warning(sprintf("No regions found for %s in run %d", gene_symbol, dnascreen_run))
    return(invisible(NULL))
  }
  
  sample_ids_region <- colnames(filtered_counts)[6:ncol(filtered_counts)]
  count_matrix_region <- filtered_counts[, sample_ids_region] %>%
    mutate_all(as.numeric)
  
  calculate_avg_percent_diff <- function(sample_col, reference_matrix) {
    reference_median <- apply(reference_matrix, 1, median, na.rm = TRUE)
    percent_diff <- abs(sample_col - reference_median) / reference_median * 100
    mean(percent_diff, na.rm = TRUE)
  }
  
  normalized_matrix <- my.count.dafr[, sample_ids_region]
  
  avg_percent_diff_vec <- c()
  for (sample in sample_ids_region) {
    reference_matrix <- normalized_matrix[, setdiff(sample_ids_region, sample), drop = FALSE]
    avg_percent_diff_vec[sample] <- calculate_avg_percent_diff(normalized_matrix[[sample]], reference_matrix)
  }
  
  avg_percent_diff_vec_filtered <- avg_percent_diff_vec[names(avg_percent_diff_vec) != target_sample_dots]
  top20_samples <- names(sort(avg_percent_diff_vec_filtered))[1:min(20, length(avg_percent_diff_vec_filtered))]
  
  target_vector <- count_matrix_region[[target_sample_dots]]
  region_medians_top20 <- apply(count_matrix_region[, top20_samples, drop = FALSE],
                                1, median, na.rm = TRUE)
  read_ratio_vector <- round(target_vector / region_medians_top20, 2)
  
  plot_df <- filtered_counts %>%
    mutate(read_ratio = read_ratio_vector) %>%
    dplyr::select(start, end, read_ratio, exon) %>%
    dplyr::arrange(start)
  
  step_df <- data.frame()
  for (i in 1:(nrow(plot_df) - 1)) {
    curr <- plot_df[i, ]
    next_row <- plot_df[i + 1, ]
    mid <- (curr$end + next_row$start) / 2
    
    step_df <- rbind(step_df,
                     data.frame(x = curr$start, y = curr$read_ratio),
                     data.frame(x = mid, y = curr$read_ratio),
                     data.frame(x = mid, y = next_row$read_ratio))
  }
  step_df <- rbind(step_df,
                   data.frame(x = plot_df$start[nrow(plot_df)], y = plot_df$read_ratio[nrow(plot_df)]),
                   data.frame(x = plot_df$end[nrow(plot_df)], y = plot_df$read_ratio[nrow(plot_df)]))
  
  has_cnv_in_gene <- (
    as.numeric(cnv_chr) == user_chr &&
      cnv_end >= user_start &&
      cnv_start <= user_end
  )
  
  rr_plot <- ggplot(plot_df) +
    {if (has_cnv_in_gene)
      geom_rect(aes(xmin = cnv_start, xmax = cnv_end, ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE, fill = "blue", alpha = 0.01)} +
    geom_segment(aes(x = start, xend = end, y = read_ratio, yend = read_ratio),
                 color = "#1f77b4", linewidth = 1.2) +
    geom_path(data = step_df, aes(x = x, y = y), color = "#1f77b4", linewidth = 0.4) +
    geom_hline(yintercept = 1,    linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 1.5,  linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0.5,  linetype = "dashed", color = "red") +
    geom_hline(yintercept = 1.65, linetype = "dotted", color = "orange") +
    geom_hline(yintercept = 0.35, linetype = "dotted", color = "orange") +
    geom_hline(yintercept = 1.35, linetype = "dotted", color = "orange") +
    geom_hline(yintercept = 0.65, linetype = "dotted", color = "orange") +
    annotate("text", x = user_start, y = 1.5,
             label = "Read Ratio = 1.5 (Duplication)",
             hjust = 0, vjust = -0.5, color = "red", size = 3.5) +
    annotate("text", x = user_start, y = 0.5,
             label = "Read Ratio = 0.5 (Deletion)",
             hjust = 0, vjust = -0.5, color = "red", size = 3.5) +
    annotate("text", x = user_start, y = 1.65,
             label = "RR = 1.65",
             hjust = 0, vjust = -0.5, color = "orange", size = 3.5) +
    annotate("text", x = user_start, y = 0.35,
             label = "RR = 0.35",
             hjust = 0, vjust = -0.5, color = "orange", size = 3.5) +
    annotate("text", x = user_start, y = 1.35,
             label = "RR = 1.35",
             hjust = 0, vjust = -0.5, color = "orange", size = 3.5) +
    annotate("text", x = user_start, y = 0.65,
             label = "RR = 0.65",
             hjust = 0, vjust = -0.5, color = "orange", size = 3.5) +
    theme_bw() +
    labs(
      title    = paste0("NGS Read Ratio for Sample ", str_extract(target_sample_dots, "DNS\\d+")),
      subtitle = paste0("chr", user_chr, ":", user_start, "-", user_end),
      x = NULL,
      y = "Read Ratio"
    ) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(
      breaks = sort(unique(c(seq(0, 4.5, by = 0.5), 0.35, 0.65, 1.35, 1.65))),
      limits = c(0, 4.5)
    ) +
    xlim(user_start, user_end) +
    theme(
      plot.title   = element_text(hjust = 0.5),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin  = margin(t = 0, r = 5, b = 0, l = 5)
    )
  
  exons_region <- exons_combined %>%
    filter(gene == gene_symbol,
           chromosome_name == as.character(user_chr),
           exon_chrom_start >= user_start,
           exon_chrom_end <= user_end)
  
  transcript_id <- clinically_relevant_transcripts[[gene_symbol]]
  exons_sorted <- exons_region %>% arrange(exon_chrom_start)
  strand <- unique(exons_sorted$strand)
  arrow_symbol <- ifelse(strand == -1, "<", ">")
  arrow_spacing <- 500
  arrow_y <- 0.5
  
  arrow_df <- data.frame()
  if (nrow(exons_sorted) > 1) {
    for (i in 1:(nrow(exons_sorted) - 1)) {
      intron_start <- exons_sorted$exon_chrom_end[i] + 1
      intron_end <- exons_sorted$exon_chrom_start[i + 1] - 1
      if (intron_start < intron_end) {
        arrow_positions <- seq(intron_start, intron_end, by = arrow_spacing)
        arrow_df <- rbind(
          arrow_df,
          data.frame(x = arrow_positions, y = arrow_y, label = arrow_symbol)
        )
      }
    }
  }
  
  gene_mid <- exons_region %>%
    summarize(mid = mean(range(exon_chrom_start, exon_chrom_end))) %>%
    pull(mid)
  
  exon_plot <- ggplot(exons_region) +
    geom_hline(yintercept = arrow_y, color = "black", linewidth = 0.4) +
    geom_text(data = arrow_df,
              aes(x = x, y = y, label = label),
              size = 3, color = "black", fontface = "bold") +
    geom_rect(aes(xmin = exon_chrom_start, xmax = exon_chrom_end,
                  ymin = 0.4, ymax = 0.6, fill = gene),
              color = "black", alpha = 1) +
    geom_text(data = data.frame(x = gene_mid, y = 0.7,
                                label = paste0(gene_symbol, " (", transcript_id, ")")),
              aes(x = x, y = y, label = label),
              size = 3, fontface = "bold") +
    theme_void() +
    scale_fill_brewer(palette = "Set2") +
    theme(
      legend.position = "none",
      plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
    ) +
    xlim(user_start, user_end)
  
  combined <- rr_plot / exon_plot +
    plot_layout(heights = c(5.8, 0.8))
  
  dns_id <- str_extract(target_sample, "DNS\\d+")
  file_name <- sprintf("%s_%s_%s_run%d.jpeg",
                       dns_id, gene_symbol, cnv_type, dnascreen_run)
  
  ggsave(
    filename = file.path(output_dir, file_name),
    plot = combined,
    width = 12,
    height = 7,
    units = "in",
    dpi = 150
  )
  
  message(sprintf("Saved: %s", file_name))
}


# ============================================================
# BATCH LOOP
# ============================================================

manual_calls <- manual_inspection_calls %>%
  # filter(acgh_status == "concordant") %>%
  dplyr::select(Sample, Gene, Chromosome, Start, End, CNV_Type, run)

for (i in seq_len(nrow(fp_calls))) {
  row <- fp_calls[i, ]
  message(sprintf("Plotting %d / %d: %s %s %s",
                  i, nrow(fp_calls),
                  str_extract(row$Sample, "DNS\\d+"),
                  row$Gene, row$CNV_Type))
  tryCatch(
    plot_rr(
      target_sample  = row$Sample,
      gene_symbol    = row$Gene,
      cnv_chr        = row$Chromosome,
      cnv_start      = row$Start,
      cnv_end        = row$End,
      cnv_type       = row$CNV_Type,
      dnascreen_run  = row$run,
      output_dir     = output_dir,
      exons_combined = exons_combined,
      clinically_relevant_transcripts = clinically_relevant_transcripts
    ),
    error = function(e) {
      message(sprintf("  ERROR row %d (%s %s): %s",
                      i,
                      str_extract(row$Sample, "DNS\\d+"),
                      row$Gene,
                      e$message))
    }
  )
}
