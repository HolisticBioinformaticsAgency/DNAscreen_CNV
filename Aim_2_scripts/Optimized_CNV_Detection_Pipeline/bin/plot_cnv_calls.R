#!/usr/bin/env Rscript
# ============================================================
# Per-call read ratio (+ VAF, + optional raw depth) plots
# Port of plot_rr_with_vaf.R. VAF status uses the framework's own rules,
# so the plot always agrees with the stage-4 decision.
# Output: calls/<final_status>/<run>_<ID>_<gene>_<type>_<caller>_chr<c>-<start>-<end>.jpeg
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(patchwork)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--run"),
  make_option("--stage4"),
  make_option("--stage2"),
  make_option("--counts"),
  make_option("--outliers"),
  make_option("--vcf_dir"),
  make_option("--bed"),
  make_option("--vcf_pattern",    default = "\\.sorted\\.vcf\\.gz$"),
  make_option("--id_regex",       default = "DNS\\d+"),
  make_option("--exon_track_rds", default = ""),
  make_option("--n_ref",          type = "integer", default = 20L),
  make_option("--vaf_max",        type = "double",  default = 0.95),
  make_option("--dp_min",         type = "double",  default = 30),
  make_option("--rr_min_max",     type = "double",  default = 2.5),
  make_option("--show_vaf",       default = "true"),
  make_option("--show_vaf_legend", default = "false"),
  make_option("--show_raw_depth", default = "false"),
  make_option("--plot_prelim_rejected", default = "false"),
  make_option("--base_size",      type = "double",  default = 14),
  make_option("--utils")
)))
source(opt$utils)
as_flag <- function(x) tolower(x) %in% c("true", "t", "1", "yes")
run <- opt$run

show_vaf        <- as_flag(opt$show_vaf)
show_vaf_legend <- as_flag(opt$show_vaf_legend)
show_raw_depth  <- as_flag(opt$show_raw_depth)

# ── Calls to plot ───────────────────────────────────────────────────────────
stage4 <- read_tsv_safe(opt$stage4)
to_plot <- if (nrow(stage4) > 0) stage4 %>% filter(role != "merged") else tibble()

if (as_flag(opt$plot_prelim_rejected)) {
  s2 <- read_tsv_safe(opt$stage2)
  if (nrow(s2) > 0) {
    rej <- s2 %>% filter(!(prelim_pass %in% TRUE)) %>%
      mutate(final_status = "rejected_prelim_filter", role = "single_caller")
    to_plot <- bind_rows(to_plot, rej)
  }
}
dir.create("calls", showWarnings = FALSE)
if (nrow(to_plot) == 0) {
  message(sprintf("%s: no calls to plot", run))
  quit(save = "no")
}

# ── Shared inputs ───────────────────────────────────────────────────────────
counts   <- load_counts(opt$counts)
outliers <- read_tsv_safe(opt$outliers)
retained <- intersect(outliers$sample[!(outliers$is_outlier %in% TRUE)], colnames(counts$raw))
norm     <- cpt_normalise(counts$raw[, retained, drop = FALSE])
apd      <- compute_apd(norm)
bed      <- read_framework_bed(opt$bed)
vcf_lookup <- build_vcf_lookup(opt$vcf_dir, opt$id_regex, opt$vcf_pattern)

exon_track <- if (nzchar(opt$exon_track_rds) && file.exists(opt$exon_track_rds)) {
  readRDS(opt$exon_track_rds) %>%
    transmute(gene, chr = std_chr(chromosome_name),
              start = exon_chrom_start, end = exon_chrom_end,
              strand = if_else(strand == -1, "-", "+"))
} else {
  bed %>% transmute(gene, chr, start, end,
                    strand = unname(coalesce(gene_strand[gene], "+")))
}

pt2mm       <- function(pt) pt / 2.845
base_size   <- opt$base_size
annot_major <- pt2mm(base_size * 0.85)
annot_minor <- pt2mm(base_size * 0.75)
annot_arrow <- pt2mm(base_size * 0.75)

base_theme <- theme_bw(base_size = base_size) +
  theme(
    plot.title    = element_text(hjust = 0.5, size = base_size, face = "bold"),
    plot.subtitle = element_text(hjust = 0,   size = base_size * 0.75),
    axis.title    = element_text(size = base_size * 0.9),
    axis.text     = element_text(size = base_size * 0.8),
    legend.title  = element_text(size = base_size * 0.8),
    legend.text   = element_text(size = base_size * 0.75)
  )

# ── One plot ────────────────────────────────────────────────────────────────
plot_call <- function(cnv) {
  chr       <- std_chr(cnv$Chromosome)
  cnv_start <- cnv$Start
  cnv_end   <- cnv$End
  gene      <- cnv$Gene
  dns_id    <- coalesce(str_extract(cnv$sample_std, opt$id_regex), cnv$sample_std)

  gene_ex <- exon_track %>% filter(gene == !!gene, chr == !!chr)
  if (nrow(gene_ex) == 0) {
    user_start <- cnv_start - 5000; user_end <- cnv_end + 5000
  } else {
    user_start <- min(gene_ex$start, cnv_start); user_end <- max(gene_ex$end, cnv_end)
  }

  # Read ratios for every exon in the gene
  rr <- exon_read_ratios(cnv, counts, norm, apd, n_ref = opt$n_ref)
  if (nrow(rr) == 0) {
    warning(sprintf("No gene exons in counts for %s", cnv$call_id)); return(invisible(NULL))
  }
  plot_df <- rr %>% transmute(start = exon_start, end = exon_end, read_ratio) %>% arrange(start)
  user_start <- min(user_start, plot_df$start); user_end <- max(user_end, plot_df$end)
  x_range <- user_end - user_start
  label_x <- user_start + x_range * 0.015

  max_rr      <- max(plot_df$read_ratio[is.finite(plot_df$read_ratio)], 0, na.rm = TRUE)
  rr_plot_max <- max(ceiling(max_rr / 0.5) * 0.5, opt$rr_min_max)
  plot_df$read_ratio <- pmin(plot_df$read_ratio, rr_plot_max)

  step_df <- NULL
  n <- nrow(plot_df)
  if (n > 1) {
    step_df <- bind_rows(lapply(seq_len(n - 1), function(i) {
      mid <- (plot_df$end[i] + plot_df$start[i + 1]) / 2
      tibble(x = c(plot_df$start[i], mid, mid),
             y = c(plot_df$read_ratio[i], plot_df$read_ratio[i], plot_df$read_ratio[i + 1]))
    }))
  }
  step_df <- bind_rows(step_df, tibble(x = c(plot_df$start[n], plot_df$end[n]),
                                       y = rep(plot_df$read_ratio[n], 2)))

  plot_title    <- paste0("NGS Read Ratio for BAM: ", dns_id, " (", run, ")")
  partner_txt   <- if (!is.na(cnv$partner_call_id %||% NA)) " | dual-caller" else ""
  plot_subtitle <- paste0(gene, " | chr", chr, ":", format(cnv_start, big.mark = ","), "-",
                          format(cnv_end, big.mark = ","), " | ", cnv$CNV_Type, " | ",
                          cnv$caller, partner_txt, "\nStatus: ", cnv$final_status,
                          if (!is.null(cnv$rr_reason) && !is.na(cnv$rr_reason))
                            paste0(" | RR ", cnv$rr_reason, " | VAF ", cnv$vaf_reason) else "")

  cnv_rect <- geom_rect(data = tibble(xmin = cnv_start, xmax = cnv_end),
                        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                        inherit.aes = FALSE, fill = "blue", alpha = 0.08)

  # ── Raw depth panel ──
  raw_depth_plot <- NULL
  if (show_raw_depth) {
    bam_path <- counts$bams[[cnv$sample_std]]
    raw_depth_plot <- tryCatch({
      seqn <- if (paste0("chr", chr) %in% seqnames(seqinfo(BamFile(bam_path)))) paste0("chr", chr) else chr
      roi  <- GRanges(seqn, IRanges(user_start, user_end))
      reads <- GenomicAlignments::readGAlignments(BamFile(bam_path), param = ScanBamParam(which = roi))
      cov  <- as.numeric(GenomicAlignments::coverage(reads)[[seqn]][user_start:user_end])
      ggplot(tibble(position = seq(user_start, user_end), depth = cov), aes(position, depth)) +
        geom_line(color = "#009E73", linewidth = 0.5) +
        scale_x_continuous(labels = scales::comma, limits = c(user_start, user_end)) +
        scale_y_continuous(labels = scales::comma) +
        base_theme +
        labs(title = plot_title, subtitle = plot_subtitle, x = NULL, y = "Depth") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              plot.margin = margin(5, 5, 0, 5))
    }, error = function(e) { warning(conditionMessage(e)); NULL })
  }

  # ── Read ratio panel ──
  rr_hlines <- list(
    list(y = 1.65, lty = "dotted", col = "orange", lbl = "RR = 1.65",              sz = annot_minor),
    list(y = 1.5,  lty = "dashed", col = "red",    lbl = "RR = 1.5 (Duplication)", sz = annot_major),
    list(y = 1.35, lty = "dotted", col = "orange", lbl = "RR = 1.35",              sz = annot_minor),
    list(y = 0.65, lty = "dotted", col = "orange", lbl = "RR = 0.65",              sz = annot_minor),
    list(y = 0.5,  lty = "dashed", col = "red",    lbl = "RR = 0.5 (Deletion)",    sz = annot_major),
    list(y = 0.35, lty = "dotted", col = "orange", lbl = "RR = 0.35",              sz = annot_minor)
  )
  rr_plot <- ggplot(plot_df) +
    cnv_rect +
    geom_segment(aes(x = start, xend = end, y = read_ratio, yend = read_ratio),
                 color = "#1f77b4", linewidth = 1.2) +
    geom_path(data = step_df, aes(x = x, y = y), color = "#1f77b4", linewidth = 0.4) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    lapply(rr_hlines, function(h) geom_hline(yintercept = h$y, linetype = h$lty, color = h$col)) +
    lapply(rr_hlines, function(h) annotate("text", x = label_x, y = h$y, label = h$lbl,
                                           hjust = 0, vjust = -0.4, color = h$col, size = h$sz)) +
    scale_x_continuous(labels = scales::comma, limits = c(user_start, user_end)) +
    scale_y_continuous(breaks = sort(unique(c(seq(0, rr_plot_max, by = 0.5), 0.35, 0.65, 1.35, 1.65))),
                       limits = c(0, rr_plot_max)) +
    base_theme +
    labs(title    = if (is.null(raw_depth_plot)) plot_title else NULL,
         subtitle = if (is.null(raw_depth_plot)) plot_subtitle else NULL,
         x = NULL, y = "Read Ratio") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(5, 5, 0, 5))

  # ── VAF panel (whole gene window, statuses from the call's CNV type) ──
  vaf_plot <- NULL
  if (show_vaf) {
    vcf <- find_vcf(cnv$sample_std, vcf_lookup, opt$id_regex)
    vaf_df <- extract_snvs_in_region(vcf, chr, user_start, user_end, opt$dp_min) %>%
      filter(vaf > 0, vaf < opt$vaf_max) %>%
      annotate_vaf_status(cnv$CNV_Type) %>%
      add_wilson_bounds() %>%
      mutate(cnv_location = if_else(pos >= cnv_start & pos <= cnv_end, "within call", "outside call"))

    if (nrow(vaf_df) > 0) {
      vaf_plot <- ggplot(vaf_df, aes(x = pos, y = vaf, colour = vaf_status, shape = cnv_location)) +
        cnv_rect +
        geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                      width = x_range * 0.005, linewidth = 0.4, alpha = 0.5) +
        geom_point(size = 2.5, alpha = 0.85) +
        geom_hline(yintercept = 0.5,         linetype = "dashed", color = "gray50") +
        geom_hline(yintercept = c(1/3, 2/3), linetype = "dotted", color = "gray40") +
        annotate("text", x = label_x, y = 1/3, label = "VAF = 0.33 (consistent with duplication)",
                 hjust = 0, vjust = -0.4, color = "gray40", size = annot_minor) +
        annotate("text", x = label_x, y = 2/3, label = "VAF = 0.67 (consistent with duplication)",
                 hjust = 0, vjust = -0.4, color = "gray40", size = annot_minor) +
        scale_colour_manual(values = c(concordant = "#2ca02c", discordant = "#d62728", neutral = "black"),
                            name = "VAF status", drop = FALSE) +
        scale_shape_manual(values = c("within call" = 16, "outside call" = 1), name = "Location") +
        scale_x_continuous(labels = scales::comma, limits = c(user_start, user_end)) +
        scale_y_continuous(limits = c(0, 1), breaks = c(0, 1/3, 0.5, 2/3, 1),
                           labels = c("0", "0.33", "0.50", "0.67", "1.0")) +
        base_theme +
        labs(x = "Genomic position (bp)", y = "VAF") +
        theme(legend.position = if (show_vaf_legend) "right" else "none",
              plot.margin = margin(0, 5, 0, 5))
    } else {
      msg <- if (is.na(vcf)) "No VCF found for sample" else
        sprintf("No informative SNVs (0 < VAF < %.2f, DP ≥ %d)", opt$vaf_max, as.integer(opt$dp_min))
      vaf_plot <- ggplot() +
        annotate("text", x = (user_start + user_end) / 2, y = 0.5, label = msg,
                 size = annot_major, color = "gray50") +
        scale_x_continuous(limits = c(user_start, user_end)) +
        scale_y_continuous(limits = c(0, 1)) +
        base_theme + labs(x = NULL, y = "VAF") +
        theme(plot.margin = margin(0, 5, 0, 5))
    }
  }

  # ── Gene / exon track ──
  exons_sorted <- gene_ex %>% filter(start >= user_start, end <= user_end) %>% arrange(start)
  arrow_symbol <- if (any(exons_sorted$strand == "-")) "<" else ">"
  arrow_df <- tibble()
  if (nrow(exons_sorted) > 1) {
    for (i in seq_len(nrow(exons_sorted) - 1)) {
      s <- exons_sorted$end[i] + 1; e <- exons_sorted$start[i + 1] - 1
      if (s < e) arrow_df <- bind_rows(arrow_df, tibble(x = seq(s, e, by = 500), y = 0.5, label = arrow_symbol))
    }
  }
  transcript <- unname(clinically_relevant_transcripts[gene])
  gene_label <- if (is.na(transcript)) gene else paste0(gene, " (", transcript, ")")
  exon_plot <- ggplot(exons_sorted) +
    geom_hline(yintercept = 0.5, color = "black", linewidth = 0.4) +
    {if (nrow(arrow_df) > 0) geom_text(data = arrow_df, aes(x = x, y = y, label = label),
                                       size = annot_arrow, fontface = "bold")} +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0.4, ymax = 0.6),
              fill = "#66C2A5", color = "black") +
    annotate("text", x = (user_start + user_end) / 2, y = 0.75, label = gene_label,
             size = annot_major, fontface = "bold") +
    scale_x_continuous(labels = scales::comma, limits = c(user_start, user_end)) +
    scale_y_continuous(limits = c(0.25, 1.05)) +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(2, 5, 4, 5))

  # ── Combine ──
  panels  <- list(); heights <- c()
  if (!is.null(raw_depth_plot)) { panels <- c(panels, list(raw_depth_plot)); heights <- c(heights, 2) }
  panels  <- c(panels, list(rr_plot));  heights <- c(heights, 5)
  if (!is.null(vaf_plot))       { panels <- c(panels, list(vaf_plot));       heights <- c(heights, 2.5) }
  panels  <- c(panels, list(exon_plot)); heights <- c(heights, 1.5)
  combined <- Reduce(`/`, panels) + plot_layout(heights = heights)

  out_dir <- file.path("calls", cnv$final_status)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  fname <- sprintf("%s_%s_%s_%s_%s_chr%s-%s-%s.jpeg", run, dns_id, gene, cnv$CNV_Type, cnv$caller,
                   chr, format(cnv_start, scientific = FALSE), format(cnv_end, scientific = FALSE))
  ggsave(file.path(out_dir, fname), combined, width = 12,
         height = max(4, sum(heights) / 11 * 11), units = "in", dpi = 150)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

for (i in seq_len(nrow(to_plot))) {
  cnv <- to_plot[i, ]
  tryCatch(plot_call(cnv), error = function(e)
    message(sprintf("  ERROR plotting %s: %s", cnv$call_id, conditionMessage(e))))
}
message(sprintf("%s: plotted %d calls", run, nrow(to_plot)))
