library(ExomeDepth)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(tidyr)
library(readr)
library(GenomicRanges)
library(IRanges)
library(Rsamtools)
library(GenomicAlignments)
library(VariantAnnotation)

# ============================================================
# OUTPUT DIRECTORY
# ============================================================

output_dir <- "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_2/Optimized_CNV_Detection_Pipeline/test"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# BAM DIRECTORY
# ============================================================

bam_dir <- "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/Validations/bam_and_vcf"

# ============================================================
# LOAD CACHED EXON DATA (gene track)
# ============================================================

exons_combined <- readRDS("/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_2/CNV-plot/cached_exon_data.rds")

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
# HELPER: WILSON SCORE 95% CI
# ============================================================

wilson_ci <- function(k, n, conf = 0.95) {
  z      <- qnorm(1 - (1 - conf) / 2)
  p_hat  <- k / n
  denom  <- 1 + z^2 / n
  centre <- (p_hat + z^2 / (2 * n)) / denom
  margin <- (z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))) / denom
  list(lower = centre - margin, upper = centre + margin)
}

vaf_within_ci <- function(vaf, dp, targets = c(1/3, 2/3)) {
  any(sapply(targets, function(p_null) {
    k  <- round(vaf * dp)
    ci <- wilson_ci(k, dp)
    p_null >= ci$lower & p_null <= ci$upper
  }))
}

vaf_within_ci_het <- function(vaf, dp, target = 0.5) {
  k  <- round(vaf * dp)
  ci <- wilson_ci(k, dp)
  target >= ci$lower & target <= ci$upper
}

# ============================================================
# LOCATE VCF FILES
# ============================================================

vcf_dir   <- "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/Validations/bam_and_vcf"
vcf_files <- list.files(vcf_dir, pattern = "\\.sorted\\.vcf\\.gz$", full.names = TRUE)
vcf_lookup <- setNames(
  vcf_files,
  str_extract(basename(vcf_files), "DNS\\d+")
)

# ============================================================
# EXTRACT VAFs WITHIN CNV REGION
# ============================================================

extract_vafs_in_region <- function(vcf_path, chrom, start, end) {
  
  region <- GRanges(
    seqnames = paste0("chr", gsub("^chr", "", as.character(chrom))),
    ranges   = IRanges(start = start, end = end)
  )
  
  tryCatch({
    vcf <- readVcf(TabixFile(vcf_path), genome = "hg38", param = ScanVcfParam(which = region))
    
    if (length(vcf) == 0) {
      seqlevelsStyle(region) <- "NCBI"
      vcf <- readVcf(TabixFile(vcf_path), genome = "hg38", param = ScanVcfParam(which = region))
    }
    
    if (length(vcf) == 0) return(NULL)
    
    type_field <- info(vcf)$TYPE
    if (!is.null(type_field)) {
      vcf <- vcf[as.character(type_field) == "SNV"]
    } else {
      is_snv <- width(ref(vcf)) == 1 &
        sapply(alt(vcf), function(a) all(width(a) == 1))
      vcf <- vcf[is_snv]
    }
    
    if (length(vcf) == 0) return(NULL)
    
    af_info <- info(vcf)$AF
    if (!is.null(af_info)) {
      vafs <- unlist(af_info)
    } else {
      ad         <- geno(vcf)$AD
      if (is.null(ad)) return(NULL)
      ad_mat     <- do.call(rbind, lapply(ad[, 1], function(x) x))
      ref_counts <- ad_mat[, 1]
      alt_counts <- rowSums(ad_mat[, -1, drop = FALSE])
      dp_total   <- ref_counts + alt_counts
      vafs       <- alt_counts / dp_total
      dp_vals    <- dp_total
    }
    
    dp_field <- geno(vcf)$DP
    if (!is.null(dp_field)) {
      dp_vals <- as.numeric(dp_field[, 1])
    } else if (!exists("dp_vals")) {
      dp_vals <- rep(NA_real_, length(vafs))
    }
    
    pass_depth <- !is.na(dp_vals) & dp_vals >= 30
    if (!any(pass_depth)) return(NULL)
    vafs    <- vafs[pass_depth]
    dp_vals <- dp_vals[pass_depth]
    pos_vec <- start(rowRanges(vcf))[pass_depth]
    
    tibble(vaf = as.numeric(vafs), dp = as.numeric(dp_vals), pos = pos_vec)
    
  }, error = function(e) {
    warning(sprintf("VCF read error at %s:%d-%d — %s", chrom, start, end, e$message))
    NULL
  })
}

# ============================================================
# ANNOTATE VAF STATUS
# ============================================================

annotate_vaf_status <- function(df, cnv_type) {
  df %>%
    rowwise() %>%
    mutate(
      vaf_status = if (tolower(cnv_type) == "deletion") {
        case_when(
          vaf_within_ci_het(vaf, dp, target = 0.5) &
            !vaf_within_ci_het(vaf, dp, target = 0.0) &
            !vaf_within_ci_het(vaf, dp, target = 1.0)  ~ "discordant",
          TRUE                                          ~ "neutral"
        )
      } else {
        case_when(
          vaf_within_ci_het(vaf, dp, target = 0.5) &
            vaf_within_ci(vaf, dp, targets = c(1/3, 2/3))  ~ "neutral",
          vaf_within_ci(vaf, dp, targets = c(1/3, 2/3))    ~ "concordant",
          TRUE                                              ~ "discordant"
        )
      }
    ) %>%
    ungroup()
}

# ============================================================
# MAIN PLOTTING FUNCTION
# ============================================================

plot_rr_vaf <- function(target_sample, gene_symbol, cnv_chr,
                        cnv_start, cnv_end, cnv_type,
                        dnascreen_run, output_dir,
                        exons_combined,
                        clinically_relevant_transcripts,
                        vcf_lookup,
                        bam_dir        = NULL,
                        show_raw_depth = TRUE,
                        show_rr        = TRUE,
                        show_vaf       = TRUE,
                        show_vaf_legend = TRUE,   # <-- new: toggle VAF panel legend
                        rr_plot_max    = NULL,
                        rr_min_max     = 2.0,
                        base_size      = 14) {
  
  # ── Font size helpers ────────────────────────────────────────────────────────
  pt2mm       <- function(pt) pt / 2.845
  annot_major <- pt2mm(base_size * 0.85)
  annot_minor <- pt2mm(base_size * 0.75)
  annot_arrow <- pt2mm(base_size * 0.75)
  
  # ── Load read counts ────────────────────────────────────────────────────────
  read_count_dir   <- "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/exome_depth/RData"
  read_count_rdata <- sprintf("%s/run%d_bams_9genes_25bp.fix.sorted.RData",
                              read_count_dir, dnascreen_run)
  if (!file.exists(read_count_rdata)) {
    warning(sprintf("RData not found for run %d, skipping.", dnascreen_run))
    return(invisible(NULL))
  }
  load(read_count_rdata)
  
  my.count.dafr <- as(my.count, "data.frame")
  sample_ids    <- colnames(my.count.dafr)[6:ncol(my.count.dafr)]
  count_matrix  <- my.count.dafr[, sample_ids] %>% mutate_all(as.numeric)
  
  sample_totals <- colSums(count_matrix, na.rm = TRUE)
  normalized_df <- as.data.frame(sweep(count_matrix, 2, sample_totals, "/") * 1e3)
  colnames(normalized_df)     <- sample_ids
  my.count.dafr[, sample_ids] <- normalized_df
  my.count.dafr$chromosome    <- as.numeric(gsub("^chr", "", my.count.dafr$chromosome))
  
  target_sample_with_suffix <- paste0(target_sample, ".hq.sorted.marked.bam")
  target_sample_dots        <- gsub("-", ".", target_sample_with_suffix)
  
  if (!target_sample_dots %in% sample_ids) {
    warning(sprintf("Sample %s not found in run %d, skipping.", target_sample_dots, dnascreen_run))
    return(invisible(NULL))
  }
  
  # ── Genomic window ──────────────────────────────────────────────────────────
  gene_exon_coords <- exons_combined %>%
    filter(gene == gene_symbol, chromosome_name == as.character(cnv_chr))
  
  if (nrow(gene_exon_coords) == 0) {
    user_start <- cnv_start - 5000
    user_end   <- cnv_end   + 5000
  } else {
    user_start <- min(gene_exon_coords$exon_chrom_start)
    user_end   <- max(gene_exon_coords$exon_chrom_end)
  }
  user_chr <- as.numeric(cnv_chr)
  
  filtered_counts <- my.count.dafr %>%
    dplyr::filter(chromosome == user_chr, start >= user_start, end <= user_end)
  
  if (nrow(filtered_counts) == 0) {
    warning(sprintf("No regions found for %s in run %d", gene_symbol, dnascreen_run))
    return(invisible(NULL))
  }
  
  # ── Read ratio ──────────────────────────────────────────────────────────────
  sample_ids_region   <- colnames(filtered_counts)[6:ncol(filtered_counts)]
  count_matrix_region <- filtered_counts[, sample_ids_region] %>% mutate_all(as.numeric)
  normalized_matrix   <- my.count.dafr[, sample_ids_region]
  
  avg_pct_diff <- sapply(sample_ids_region, function(s) {
    ref_mat <- normalized_matrix[, setdiff(sample_ids_region, s), drop = FALSE]
    ref_med <- apply(ref_mat, 1, median, na.rm = TRUE)
    mean(abs(normalized_matrix[[s]] - ref_med) / ref_med * 100, na.rm = TRUE)
  })
  
  top20_samples <- names(sort(avg_pct_diff[names(avg_pct_diff) != target_sample_dots]))[
    1:min(20, sum(names(avg_pct_diff) != target_sample_dots))]
  target_vector <- count_matrix_region[[target_sample_dots]]
  region_med    <- apply(count_matrix_region[, top20_samples, drop = FALSE], 1, median, na.rm = TRUE)
  rr_vector     <- round(target_vector / region_med, 2)
  
  # ── Auto-compute rr_plot_max ─────────────────────────────────────────────────
  if (is.null(rr_plot_max)) {
    max_rr      <- max(rr_vector, na.rm = TRUE)
    rr_plot_max <- max(ceiling(max_rr / 0.5) * 0.5, rr_min_max)
    message(sprintf("  Auto rr_plot_max: %.2f → %.1f (floor = %.1f)", max_rr, rr_plot_max, rr_min_max))
  }
  
  plot_df <- filtered_counts %>%
    mutate(read_ratio = rr_vector) %>%
    dplyr::select(start, end, read_ratio) %>%
    dplyr::arrange(start)
  
  step_df <- data.frame()
  for (i in 1:(nrow(plot_df) - 1)) {
    mid     <- (plot_df$end[i] + plot_df$start[i + 1]) / 2
    step_df <- rbind(step_df,
                     data.frame(x = plot_df$start[i],   y = plot_df$read_ratio[i]),
                     data.frame(x = mid,                 y = plot_df$read_ratio[i]),
                     data.frame(x = mid,                 y = plot_df$read_ratio[i + 1]))
  }
  step_df <- rbind(step_df,
                   data.frame(x = plot_df$start[nrow(plot_df)], y = plot_df$read_ratio[nrow(plot_df)]),
                   data.frame(x = plot_df$end[nrow(plot_df)],   y = plot_df$read_ratio[nrow(plot_df)]))
  
  # ── VAF extraction & annotation ─────────────────────────────────────────────
  dns_id   <- str_extract(target_sample, "DNS\\d+")
  vcf_path <- vcf_lookup[dns_id]
  vaf_df   <- NULL
  raw_vafs <- NULL
  
  if (!is.na(vcf_path) && file.exists(vcf_path)) {
    raw_vafs <- extract_vafs_in_region(vcf_path,
                                       chrom = cnv_chr,
                                       start = user_start,
                                       end   = user_end)
  }
  
  if (!is.null(raw_vafs) && nrow(raw_vafs) > 0) {
    vaf_df <- raw_vafs %>%
      filter(!is.na(vaf), !is.na(dp), vaf > 0, vaf < 0.95) %>%
      annotate_vaf_status(cnv_type = tolower(cnv_type)) %>%
      mutate(
        cnv_location = if_else(pos >= cnv_start & pos <= cnv_end, "within call", "outside call"),
        ci_lower = mapply(function(v, d) { ci <- wilson_ci(round(v * d), d); ci$lower }, vaf, dp),
        ci_upper = mapply(function(v, d) { ci <- wilson_ci(round(v * d), d); ci$upper }, vaf, dp)
      )
  }
  
  has_cnv <- cnv_end >= user_start & cnv_start <= user_end
  
  # ── Plot title & subtitle ────────────────────────────────────────────────────
  plot_title    <- paste0("NGS Read Ratio for BAM: ", dns_id)
  plot_subtitle <- paste0(gene_symbol, " | chr", user_chr, ":", user_start, "-", user_end,
                          " | ", cnv_type)
  
  # ── Shared base theme ────────────────────────────────────────────────────────
  base_theme <- theme_bw(base_size = base_size) +
    theme(
      plot.title    = element_text(hjust = 0.5, size = base_size,        face = "bold"),
      plot.subtitle = element_text(hjust = 0,   size = base_size * 0.85),
      axis.title    = element_text(size = base_size * 0.9),
      axis.text     = element_text(size = base_size * 0.8),
      legend.title  = element_text(size = base_size * 0.8),
      legend.text   = element_text(size = base_size * 0.75)
    )
  
  # ── Label x-position ────────────────────────────────────────────────────────
  x_range <- user_end - user_start
  label_x <- user_start + x_range * 0.015
  
  # ── Raw depth panel ──────────────────────────────────────────────────────────
  raw_depth_plot <- NULL
  
  if (show_raw_depth && !is.null(bam_dir)) {
    bam_path  <- file.path(bam_dir, target_sample_with_suffix)
    bam_index <- paste0(bam_path, ".bai")
    
    if (file.exists(bam_path) && file.exists(bam_index)) {
      tryCatch({
        roi      <- GRanges(seqnames = paste0("chr", user_chr),
                            ranges   = IRanges(start = user_start, end = user_end))
        bam_file <- BamFile(bam_path, index = bam_index)
        reads    <- readGAlignments(bam_file, param = ScanBamParam(which = roi))
        cov      <- coverage(reads)
        cov_vec  <- cov[[paste0("chr", user_chr)]][user_start:user_end]
        
        raw_depth_df <- data.frame(
          position = seq(user_start, user_end),
          depth    = as.numeric(cov_vec)
        )
        
        raw_depth_plot <- ggplot(raw_depth_df, aes(x = position, y = depth)) +
          geom_line(color = "#009E73", linewidth = 0.5) +
          geom_hline(yintercept = 0, color = "gray80", linewidth = 0.3) +
          scale_x_continuous(labels = scales::comma, limits = c(user_start, user_end)) +
          scale_y_continuous(labels = scales::comma) +
          base_theme +
          labs(title = plot_title, subtitle = plot_subtitle, x = NULL, y = "Depth") +
          theme(
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin  = margin(t = 5, r = 5, b = 0, l = 5)
          )
        
      }, error = function(e) {
        warning(sprintf("Raw depth failed for %s — %s", target_sample, e$message))
        raw_depth_plot <<- NULL
      })
    } else {
      warning(sprintf("BAM or index not found for %s, skipping raw depth panel.", target_sample))
    }
  }
  
  # ── Read ratio panel ─────────────────────────────────────────────────────────
  rr_plot <- NULL
  
  if (show_rr) {
    rr_title    <- if (is.null(raw_depth_plot)) plot_title    else NULL
    rr_subtitle <- if (is.null(raw_depth_plot)) plot_subtitle else NULL
    
    rr_hlines <- list(
      list(y = 1.65, lty = "dotted", col = "orange", lbl = "RR = 1.65",       sz = annot_minor),
      list(y = 1.5,  lty = "dashed", col = "red",    lbl = "RR = 1.5 (Duplication)",  sz = annot_major),
      list(y = 1.35, lty = "dotted", col = "orange", lbl = "RR = 1.35",       sz = annot_minor),
      list(y = 0.65, lty = "dotted", col = "orange", lbl = "RR = 0.65",       sz = annot_minor),
      list(y = 0.5,  lty = "dashed", col = "red",    lbl = "RR = 0.5 (Deletion)",  sz = annot_major),
      list(y = 0.35, lty = "dotted", col = "orange", lbl = "RR = 0.35",       sz = annot_minor)
    )
    
    rr_plot <- ggplot(plot_df) +
      {if (has_cnv)
        geom_rect(aes(xmin = cnv_start, xmax = cnv_end, ymin = -Inf, ymax = Inf),
                  inherit.aes = FALSE, fill = "blue", alpha = 0.01)} +
      geom_segment(aes(x = start, xend = end, y = read_ratio, yend = read_ratio),
                   color = "#1f77b4", linewidth = 1.2) +
      geom_path(data = step_df, aes(x = x, y = y), color = "#1f77b4", linewidth = 0.4) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      {lapply(rr_hlines, function(h) {
        if (h$y <= rr_plot_max)
          geom_hline(yintercept = h$y, linetype = h$lty, color = h$col)
      })} +
      {lapply(rr_hlines, function(h) {
        if (h$y <= rr_plot_max)
          annotate("text", x = label_x, y = h$y, label = h$lbl,
                   hjust = 0, vjust = -0.4, color = h$col, size = h$sz)
      })} +
      scale_x_continuous(labels = scales::comma, limits = c(user_start, user_end)) +
      scale_y_continuous(
        breaks = sort(unique(c(seq(0, rr_plot_max, by = 0.5), 0.35, 0.65, 1.35, 1.65))),
        limits = c(0, rr_plot_max)
      ) +
      base_theme +
      labs(title = rr_title, subtitle = rr_subtitle, x = NULL, y = "Read Ratio") +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(t = 5, r = 5, b = 0, l = 5)
      )
  }
  
  # ── VAF panel ────────────────────────────────────────────────────────────────
  vaf_plot <- NULL
  
  if (show_vaf) {
    vaf_title    <- if (is.null(raw_depth_plot) && is.null(rr_plot)) plot_title    else NULL
    vaf_subtitle <- if (is.null(raw_depth_plot) && is.null(rr_plot)) plot_subtitle else NULL
    
    vaf_colours <- c("concordant" = "#2ca02c", "discordant" = "#d62728", "neutral" = "black")
    vaf_shapes  <- c("within call" = 16, "outside call" = 1)
    
    # ── Resolve legend position from show_vaf_legend ─────────────────────────
    vaf_legend_position <- if (isTRUE(show_vaf_legend)) "right" else "none"
    
    if (!is.null(vaf_df) && nrow(vaf_df) > 0) {
      vaf_plot <- ggplot(vaf_df, aes(x = pos, y = vaf, colour = vaf_status, shape = cnv_location)) +
        {if (has_cnv)
          geom_rect(aes(xmin = cnv_start, xmax = cnv_end, ymin = -Inf, ymax = Inf),
                    inherit.aes = FALSE, fill = "blue", alpha = 0.01)} +
        geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                      width = x_range * 0.005, linewidth = 0.4, alpha = 0.5) +
        geom_point(size = 2.5, alpha = 0.85) +
        geom_hline(yintercept = 0.5,         linetype = "dashed", color = "gray50") +
        geom_hline(yintercept = c(1/3, 2/3), linetype = "dotted", color = "gray40") +
        annotate("text", x = label_x, y = 1/3, label = "VAF = 0.33 (consistent with duplication)",
                 hjust = 0, vjust = -0.4, color = "gray40", size = annot_minor) +
        annotate("text", x = label_x, y = 2/3, label = "VAF = 0.67 (consistent with duplication)",
                 hjust = 0, vjust = -0.4, color = "gray40", size = annot_minor) +
        scale_colour_manual(values = vaf_colours, name = "VAF status", drop = FALSE) +
        scale_shape_manual(values  = vaf_shapes,  name = "Location") +
        scale_x_continuous(labels = scales::comma, limits = c(user_start, user_end)) +
        scale_y_continuous(limits = c(0, 1), breaks = c(0, 1/3, 0.5, 2/3, 1),
                           labels = c("0", "0.33", "0.50", "0.67", "1.0")) +
        base_theme +
        labs(title = vaf_title, subtitle = vaf_subtitle, x = "Genomic position (bp)", y = "VAF") +
        theme(
          legend.position = vaf_legend_position,   # <-- controlled by show_vaf_legend
          plot.margin     = margin(t = 0, r = 5, b = 0, l = 5)
        )
    } else {
      vaf_plot <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "No VAFs available (DP \u2265 30)",
                 size = annot_major, color = "gray50") +
        scale_x_continuous(limits = c(user_start, user_end)) +
        scale_y_continuous(limits = c(0, 1)) +
        base_theme +
        labs(title = vaf_title, subtitle = vaf_subtitle, x = NULL, y = "VAF") +
        theme(plot.margin = margin(t = 0, r = 5, b = 0, l = 5))
    }
  }
  
  # ── Gene/exon track ──────────────────────────────────────────────────────────
  exons_region <- exons_combined %>%
    filter(gene == gene_symbol,
           chromosome_name == as.character(user_chr),
           exon_chrom_start >= user_start,
           exon_chrom_end   <= user_end)
  
  transcript_id <- clinically_relevant_transcripts[[gene_symbol]]
  exons_sorted  <- exons_region %>% arrange(exon_chrom_start)
  strand        <- unique(exons_sorted$strand)
  arrow_symbol  <- ifelse(strand == -1, "<", ">")
  arrow_y       <- 0.5
  arrow_df      <- data.frame()
  
  if (nrow(exons_sorted) > 1) {
    for (i in 1:(nrow(exons_sorted) - 1)) {
      intron_s <- exons_sorted$exon_chrom_end[i]       + 1
      intron_e <- exons_sorted$exon_chrom_start[i + 1] - 1
      if (intron_s < intron_e) {
        arrow_df <- rbind(arrow_df,
                          data.frame(x     = seq(intron_s, intron_e, by = 500),
                                     y     = arrow_y,
                                     label = arrow_symbol))
      }
    }
  }
  
  gene_mid  <- mean(range(exons_region$exon_chrom_start, exons_region$exon_chrom_end))
  exon_plot <- ggplot(exons_region) +
    geom_hline(yintercept = arrow_y, color = "black", linewidth = 0.4) +
    {if (nrow(arrow_df) > 0)
      geom_text(data = arrow_df, aes(x = x, y = y, label = label),
                size = annot_arrow, color = "black", fontface = "bold")} +
    geom_rect(aes(xmin = exon_chrom_start, xmax = exon_chrom_end,
                  ymin = 0.4, ymax = 0.6, fill = gene),
              color = "black", alpha = 1) +
    annotate("text", x = gene_mid, y = 0.75,
             label = paste0(gene_symbol, " (", transcript_id, ")"),
             size = annot_major, fontface = "bold") +
    scale_fill_brewer(palette = "Set2") +
    scale_x_continuous(labels = scales::comma, limits = c(user_start, user_end)) +
    scale_y_continuous(limits = c(0.25, 1.05)) +
    labs(x = NULL) +
    theme_void() +
    theme(
      legend.position = "none",
      axis.title.x    = element_text(size = base_size * 0.8, margin = margin(t = 4)),
      plot.margin     = margin(t = 2, r = 5, b = 4, l = 5)
    )
  
  # ── Combine panels dynamically ───────────────────────────────────────────────
  panel_list  <- list()
  height_list <- c()
  
  if (!is.null(raw_depth_plot)) { panel_list <- c(panel_list, list(raw_depth_plot)); height_list <- c(height_list, 2)   }
  if (!is.null(rr_plot))        { panel_list <- c(panel_list, list(rr_plot));        height_list <- c(height_list, 5)   }
  if (!is.null(vaf_plot))       { panel_list <- c(panel_list, list(vaf_plot));       height_list <- c(height_list, 2.5) }
  panel_list  <- c(panel_list, list(exon_plot))
  height_list <- c(height_list, 1.5)
  
  combined <- Reduce(`/`, panel_list) + plot_layout(heights = height_list)
  
  # ── Save ─────────────────────────────────────────────────────────────────────
  full_height  <- sum(c(2, 5, 2.5, 1.5))
  total_height <- sum(height_list) / full_height * 11
  file_name    <- sprintf("%s_%s_%s_run%d.jpeg", dns_id, gene_symbol, cnv_type, dnascreen_run)
  
  ggsave(
    filename = file.path(output_dir, file_name),
    plot     = combined,
    width    = 12,
    height   = max(4, total_height),
    units    = "in", dpi = 150
  )
  message(sprintf("Saved: %s", file_name))
}

# ============================================================
# BATCH LOOP — DECoN FP calls
# ============================================================

## For fp_no_exon_near_expected

for (i in seq_len(nrow(fp_no_exon_near_expected))) {
  row <- fp_no_exon_near_expected[i, ]
  message(sprintf("Plotting %d / %d: %s %s %s",
                  i, nrow(fp_no_exon_near_expected),
                  str_extract(row$Sample, "DNS\\d+"),
                  row$Gene, row$CNV_Type))
  tryCatch(
    plot_rr_vaf(
      target_sample   = row$Sample,
      gene_symbol     = row$Gene,
      cnv_chr         = row$Chromosome,
      cnv_start       = row$CNV_Start,
      cnv_end         = row$CNV_End,
      cnv_type        = row$CNV_Type,
      dnascreen_run   = row$run,
      output_dir      = output_dir,
      exons_combined  = exons_combined,
      clinically_relevant_transcripts = clinically_relevant_transcripts,
      vcf_lookup      = vcf_lookup,
      bam_dir         = bam_dir,
      show_raw_depth  = FALSE,
      show_rr         = TRUE,
      show_vaf        = FALSE,
      show_vaf_legend = FALSE,
      rr_plot_max     = NULL,  # auto-computed per sample
      rr_min_max      = 2.0,   # y-axis never shorter than this
      base_size       = 14
    ),
    error = function(e) {
      message(sprintf("  ERROR row %d (%s %s): %s",
                      i, str_extract(row$Sample, "DNS\\d+"),
                      row$Gene, e$message))
    }
  )
}

## For fp_intragenic

for (i in seq_len(nrow(fp_intragenic))) {
  row <- fp_intragenic[i, ]
  message(sprintf("Plotting %d / %d: %s %s %s",
                  i, nrow(fp_intragenic),
                  str_extract(row$Sample, "DNS\\d+"),
                  row$Gene, row$CNV_Type))
  tryCatch(
    plot_rr_vaf(
      target_sample   = row$Sample,
      gene_symbol     = row$Gene,
      cnv_chr         = row$Chromosome,
      cnv_start       = row$CNV_Start,
      cnv_end         = row$CNV_End,
      cnv_type        = row$CNV_Type,
      dnascreen_run   = row$run,
      output_dir      = output_dir,
      exons_combined  = exons_combined,
      clinically_relevant_transcripts = clinically_relevant_transcripts,
      vcf_lookup      = vcf_lookup,
      bam_dir         = bam_dir,
      show_raw_depth  = FALSE,
      show_rr         = TRUE,
      show_vaf        = FALSE,
      show_vaf_legend = FALSE,
      rr_plot_max     = NULL,  # auto-computed per sample
      rr_min_max      = 2.0,   # y-axis never shorter than this
      base_size       = 14
    ),
    error = function(e) {
      message(sprintf("  ERROR row %d (%s %s): %s",
                      i, str_extract(row$Sample, "DNS\\d+"),
                      row$Gene, e$message))
    }
  )
}


## For fp_high_rr_spread

for (i in seq_len(nrow(fp_high_rr_spread))) {
  row <- fp_high_rr_spread[i, ]
  message(sprintf("Plotting %d / %d: %s %s %s",
                  i, nrow(fp_high_rr_spread),
                  str_extract(row$Sample, "DNS\\d+"),
                  row$Gene, row$CNV_Type))
  tryCatch(
    plot_rr_vaf(
      target_sample   = row$Sample,
      gene_symbol     = row$Gene,
      cnv_chr         = row$Chromosome,
      cnv_start       = row$CNV_Start,
      cnv_end         = row$CNV_End,
      cnv_type        = row$CNV_Type,
      dnascreen_run   = row$run,
      output_dir      = output_dir,
      exons_combined  = exons_combined,
      clinically_relevant_transcripts = clinically_relevant_transcripts,
      vcf_lookup      = vcf_lookup,
      bam_dir         = bam_dir,
      show_raw_depth  = FALSE,
      show_rr         = TRUE,
      show_vaf        = FALSE,
      show_vaf_legend = FALSE,
      rr_plot_max     = NULL,  # auto-computed per sample
      rr_min_max      = 2.0,   # y-axis never shorter than this
      base_size       = 14
    ),
    error = function(e) {
      message(sprintf("  ERROR row %d (%s %s): %s",
                      i, str_extract(row$Sample, "DNS\\d+"),
                      row$Gene, e$message))
    }
  )
}



# ============================================================
# BATCH LOOP — manual_inspection_calls
# ============================================================

for (i in seq_len(nrow(manual_inspection_calls))) {
  row <- manual_inspection_calls[i, ]
  message(sprintf("Plotting %d / %d: %s %s %s",
                  i, nrow(manual_inspection_calls),
                  str_extract(row$Sample, "DNS\\d+"),
                  row$Gene, row$CNV_Type))
  tryCatch(
    plot_rr_vaf(
      target_sample   = row$Sample,
      gene_symbol     = row$Gene,
      cnv_chr         = row$Chromosome,
      cnv_start       = row$Start,
      cnv_end         = row$End,
      cnv_type        = row$CNV_Type,
      dnascreen_run   = row$run,
      output_dir      = output_dir,
      exons_combined  = exons_combined,
      clinically_relevant_transcripts = clinically_relevant_transcripts,
      vcf_lookup      = vcf_lookup,
      bam_dir         = bam_dir,
      show_raw_depth  = FALSE,
      show_rr         = TRUE,
      show_vaf        = TRUE,
      show_vaf_legend = FALSE,
      rr_plot_max     = NULL,  # auto-computed per sample
      rr_min_max      = 2.5,   # y-axis never shorter than this
      base_size       = 14
    ),
    error = function(e) {
      message(sprintf("  ERROR row %d (%s %s): %s",
                      i, str_extract(row$Sample, "DNS\\d+"),
                      row$Gene, e$message))
    }
  )
}

# ============================================================
# BATCH LOOP — high_confidence_calls
# ============================================================

for (i in seq_len(nrow(high_confidence_calls))) {
  row <- high_confidence_calls[i, ]
  message(sprintf("Plotting %d / %d: %s %s %s",
                  i, nrow(high_confidence_calls),
                  str_extract(row$Sample, "DNS\\d+"),
                  row$Gene, row$CNV_Type))
  tryCatch(
    plot_rr_vaf(
      target_sample   = row$Sample,
      gene_symbol     = row$Gene,
      cnv_chr         = row$Chromosome,
      cnv_start       = row$Start,
      cnv_end         = row$End,
      cnv_type        = row$CNV_Type,
      dnascreen_run   = row$run,
      output_dir      = output_dir,
      exons_combined  = exons_combined,
      clinically_relevant_transcripts = clinically_relevant_transcripts,
      vcf_lookup      = vcf_lookup,
      bam_dir         = bam_dir,
      show_raw_depth  = FALSE,
      show_rr         = TRUE,
      show_vaf        = TRUE,
      show_vaf_legend = FALSE,
      rr_plot_max     = NULL,
      rr_min_max      = 2.5,
      base_size       = 14
    ),
    error = function(e) {
      message(sprintf("  ERROR row %d (%s %s): %s",
                      i, str_extract(row$Sample, "DNS\\d+"),
                      row$Gene, e$message))
    }
  )
}

# ============================================================
# BATCH LOOP — decon_unconfirmed
# ============================================================

# decon_unconfirmed <- decon_stage1_rc %>% filter(acgh_status == "tbc")
#
# for (i in seq_len(nrow(decon_unconfirmed))) {
#   row <- decon_unconfirmed[i, ]
#   message(sprintf("Plotting %d / %d: %s %s %s",
#                   i, nrow(decon_unconfirmed),
#                   str_extract(row$Sample, "DNS\\d+"),
#                   row$Gene, row$CNV_Type))
#   tryCatch(
#     plot_rr_vaf(
#       target_sample   = row$Sample,
#       gene_symbol     = row$Gene,
#       cnv_chr         = row$Chromosome,
#       cnv_start       = row$Start,
#       cnv_end         = row$End,
#       cnv_type        = row$CNV_Type,
#       dnascreen_run   = row$run,
#       output_dir      = output_dir,
#       exons_combined  = exons_combined,
#       clinically_relevant_transcripts = clinically_relevant_transcripts,
#       vcf_lookup      = vcf_lookup,
#       bam_dir         = bam_dir,
#       show_raw_depth  = TRUE,
#       show_rr         = TRUE,
#       show_vaf        = TRUE,
#       rr_plot_max     = NULL,
#       rr_min_max      = 2.0,
#       base_size       = 14
#     ),
#     error = function(e) {
#       message(sprintf("  ERROR row %d (%s %s): %s",
#                       i, str_extract(row$Sample, "DNS\\d+"),
#                       row$Gene, e$message))
#     }
#   )
# }