library(dplyr)
library(GenomicRanges)
library(stringr)
library(ExomeDepth)
library(tidyr)
library(ggplot2)
library(VariantAnnotation)

# ============================================================
# STANDARDISE FUNCTIONS
# ============================================================

standardise_decon <- function(df) {
  df %>%
    transmute(
      Sample,
      sample_id_trimmed = sub(".*-(DNS\\d+)_S\\d+.*", "\\1", Sample),
      run = as.numeric(gsub("run", "", run)),
      Chromosome = as.character(Chromosome),
      Start = Start,
      End = End,
      Gene = Gene,
      CNV_Type = CNV.type,
      Quality = BF
    )
}

standardise_clearcnv <- function(df) {
  df %>%
    transmute(
      Sample = sample,
      sample_id_trimmed = sub(".*-(DNS\\d+)_S\\d+.*", "\\1", sample),
      run = as.numeric(gsub("run", "", run)),
      Chromosome = as.character(chr),
      Start = start,
      End = end,
      Gene = gene,
      CNV_Type = recode(aberration, "DEL" = "deletion", "DUP" = "duplication"),
      Quality = score
    )
}

# ============================================================
# LOAD FILES
# ============================================================

decon_file <- read.csv(
  "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_2/Optimized_CNV_Detection_Pipeline/output/combined_cnv_calls_decon_exon_bed.tsv",
  sep = "\t"
)

clearcnv_file <- read.csv(
  "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_2/Optimized_CNV_Detection_Pipeline/output/combined_cnv_calls_clearcnv_exon_bed_001.tsv",
  sep = "\t"
)

all_samples <- read.csv(
  "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/CI_plots/all_samples_combined.csv"
)

# ============================================================
# STANDARDISE CALLS
# ============================================================

perrun_decon_calls <- standardise_decon(decon_file)
perrun_decon_calls$Sample <- sub("\\.hq\\.sorted\\.marked.*", "", perrun_decon_calls$Sample)

perrun_clearcnv_calls <- standardise_clearcnv(clearcnv_file)
perrun_clearcnv_calls$Sample <- sub("\\.hq\\.sorted\\.marked.*", "", perrun_clearcnv_calls$Sample)

# ============================================================
# ASSIGN aCGH STATUS
# (assign_acgh_status, assign_non_vscnv, assign_cyto_overlap
#  are assumed to be already defined/loaded in your environment)
# ============================================================

decon_with_acgh <- assign_acgh_status(perrun_decon_calls, acgh_ran_samples)
decon_with_acgh <- assign_cyto_overlap(decon_with_acgh, cyto_master_table)

clearcnv_with_acgh <- assign_acgh_status(perrun_clearcnv_calls, acgh_ran_samples)
clearcnv_with_acgh <- assign_cyto_overlap(clearcnv_with_acgh, cyto_master_table)

# ============================================================
# COUNTING CALLS WITH CERTAIN FEATURES (RR AND VAF)
# ============================================================

get_decon_exon_read_ratios <- function(decon_calls) {
  
  runs_needed <- unique(decon_calls$run)
  result_list <- list()
  
  for (run_num in runs_needed) {
    
    rdata_path <- sprintf(
      "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/exome_depth/RData/run%d_bams_9genes_25bp.fix.sorted.RData",
      run_num
    )
    
    if (!file.exists(rdata_path)) {
      warning(sprintf("Run %d: RData not found, skipping.", run_num))
      next
    }
    
    load(rdata_path)  # loads my.count
    
    my.count.dafr <- as(my.count, "data.frame")
    sample_ids    <- colnames(my.count.dafr)[6:ncol(my.count.dafr)]
    count_matrix  <- my.count.dafr[, sample_ids] %>%
      mutate_all(as.numeric)
    
    # CPT normalisation only (no GC loess)
    sample_totals     <- colSums(count_matrix, na.rm = TRUE)
    normalized_counts <- sweep(count_matrix, 2, sample_totals, FUN = "/") * 1e3
    normalized_df     <- as.data.frame(normalized_counts)
    colnames(normalized_df) <- sample_ids
    
    my.count.dafr[, sample_ids] <- normalized_df
    my.count.dafr$chromosome    <- gsub("^chr", "", my.count.dafr$chromosome)
    my.count.dafr$chromosome    <- as.numeric(my.count.dafr$chromosome)
    
    # Gene annotation from exon column e.g. "PCSK9_cds_0" -> "PCSK9"
    my.count.dafr$Gene <- str_extract(as.character(my.count.dafr$exon), "^[^_]+")
    
    # name_map: original col name -> normalised sample name (dashes, no suffix)
    sample_ids_normalised <- sample_ids %>%
      str_remove("\\.hq\\.sorted\\.marked.*$") %>%
      str_replace_all("\\.", "-")
    name_map <- setNames(sample_ids_normalised, sample_ids)
    
    # APD-based top 20 reference selection (computed once per run, across all exons)
    avg_percent_diff_vec <- sapply(sample_ids, function(s) {
      ref        <- my.count.dafr[, setdiff(sample_ids, s), drop = FALSE]
      ref_median <- apply(ref, 1, median, na.rm = TRUE)
      mean(abs(my.count.dafr[[s]] - ref_median) /
             (ref_median + 1e-9) * 100, na.rm = TRUE)
    })
    names(avg_percent_diff_vec) <- sample_ids
    
    run_calls <- decon_calls %>% filter(run == run_num)
    
    run_results <- lapply(seq_len(nrow(run_calls)), function(i) {
      cnv             <- run_calls[i, ]
      cnv_chr         <- as.character(cnv$Chromosome)
      cnv_start       <- cnv$Start
      cnv_end         <- cnv$End
      cnv_gene        <- cnv$Gene
      cnv_sample_norm <- cnv$Sample
      
      target_col <- names(name_map)[name_map == cnv_sample_norm]
      if (length(target_col) != 1) {
        warning(sprintf("Sample %s not found in run %d", cnv_sample_norm, run_num))
        return(NULL)
      }
      
      # Top 20 references excluding target
      candidates    <- setdiff(sample_ids, target_col)
      top20_samples <- names(sort(avg_percent_diff_vec[candidates]))[
        1:min(20, length(candidates))]
      
      # All exons in this gene and chromosome
      gene_mask  <- my.count.dafr$Gene       == cnv_gene &
        my.count.dafr$chromosome == as.numeric(cnv_chr)
      gene_exons <- my.count.dafr[gene_mask, ] %>%
        arrange(start) %>%
        mutate(row_in_gene = row_number())
      
      if (nrow(gene_exons) == 0) return(NULL)
      
      # Tag each exon relative to the CNV call
      in_call_mask <- gene_exons$start <= cnv_end &
        gene_exons$end   >= cnv_start
      
      if (!any(in_call_mask)) return(NULL)
      
      in_idx  <- which(in_call_mask)
      min_idx <- min(in_idx)
      max_idx <- max(in_idx)
      
      exon_type <- case_when(
        in_call_mask                             ~ "in_call",
        gene_exons$row_in_gene == (min_idx - 1) ~ "adjacent",
        gene_exons$row_in_gene == (max_idx + 1) ~ "adjacent",
        TRUE                                     ~ "other"
      )
      
      # Compute read ratios using top 20 references
      target_counts <- as.numeric(gene_exons[[target_col]])
      ref_medians   <- apply(gene_exons[, top20_samples, drop = FALSE],
                             1, median, na.rm = TRUE)
      exon_rr       <- round(target_counts / ref_medians, 2)
      
      tibble(
        Sample            = cnv$Sample,
        sample_id_trimmed = cnv$sample_id_trimmed,
        run               = cnv$run,
        Gene              = cnv_gene,
        CNV_Type          = cnv$CNV_Type,
        Chromosome        = cnv$Chromosome,
        CNV_Start         = cnv_start,
        CNV_End           = cnv_end,
        acgh_status       = cnv$acgh_status,
        exon_start        = gene_exons$start,
        exon_end          = gene_exons$end,
        exon_type         = exon_type,
        read_count        = target_counts,
        read_ratio        = exon_rr
      )
    })
    
    result_list[[as.character(run_num)]] <- bind_rows(run_results)
  }
  
  bind_rows(result_list)
}

decon_with_acgh_assessed <- decon_with_acgh %>%
  filter(!acgh_status %in% c("tbc", "inconclusive"))

# Run on all DECoN calls
decon_rr_exons <- get_decon_exon_read_ratios(decon_with_acgh_assessed)



# ============================================================
# COUNTING CALLS WITH CERTAIN RR FEATURES
# ============================================================

fp_no_exon_near_expected <- decon_rr_exons %>%
  filter(acgh_status == "discordant", exon_type == "in_call") %>%
  mutate(
    expected_rr      = if_else(CNV_Type == "deletion", 0.5, 1.5),
    near_expected    = abs(read_ratio - expected_rr) <= 0.15
  ) %>%
  group_by(Sample, run, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type) %>%
  summarise(
    n_exons              = n(),
    any_near_expected    = any(near_expected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!any_near_expected)

n_fp_no_exon_near_expected <- nrow(fp_no_exon_near_expected)
cat("FP calls where no in-call exon is within ±0.15 of expected RR:",
    n_fp_no_exon_near_expected, "\n")


# Intragenic events: any exon with RR <=0.35 AND >=1.65
fp_intragenic <- decon_rr_exons %>%
  filter(acgh_status == "discordant") %>%
  group_by(Sample, run, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type) %>%
  summarise(
    has_intragenic_event = any(read_ratio >= 1.35, na.rm = TRUE) &
      any(read_ratio <= 0.65, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(has_intragenic_event)

n_fp_intragenic <- nrow(fp_intragenic)
cat("FP calls with intragenic events (both RR >= 1.35 and RR <= 0.65):",
    n_fp_intragenic, "\n")

# High RR spread within call: max - min RR > 0.5
fp_high_rr_spread <- decon_rr_exons %>%
  filter(acgh_status == "discordant", exon_type == "in_call") %>%
  group_by(Sample, run, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type) %>%
  summarise(
    n_exons    = n(),
    max_rr     = max(read_ratio, na.rm = TRUE),
    min_rr     = min(read_ratio, na.rm = TRUE),
    rr_spread  = max_rr - min_rr,
    .groups = "drop"
  ) %>%
  filter(n_exons > 1, rr_spread > 0.5)

n_fp_high_rr_spread <- nrow(fp_high_rr_spread)
cat("FP calls where in-call exon RR spread (max - min) > 0.5:",
    n_fp_high_rr_spread, "\n")



# ============================================================
# HELPER: WILSON SCORE 95% CI
# ============================================================

wilson_ci <- function(k, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  p_hat <- k / n
  denom <- 1 + z^2 / n
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

vcf_dir <- "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/Validations/bam_and_vcf"

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
    
    # ── SNV-only filter ───────────────────────────────────────────────────────
    type_field <- info(vcf)$TYPE
    if (!is.null(type_field)) {
      vcf <- vcf[as.character(type_field) == "SNV"]
    } else {
      is_snv <- width(ref(vcf)) == 1 &
        sapply(alt(vcf), function(a) all(width(a) == 1))
      vcf <- vcf[is_snv]
    }
    
    if (length(vcf) == 0) return(NULL)
    
    # ── VAF extraction ────────────────────────────────────────────────────────
    af_info <- info(vcf)$AF
    
    if (!is.null(af_info)) {
      vafs <- unlist(af_info)
    } else {
      ad <- geno(vcf)$AD
      if (is.null(ad)) return(NULL)
      ad_mat     <- do.call(rbind, lapply(ad[, 1], function(x) x))
      ref_counts <- ad_mat[, 1]
      alt_counts <- rowSums(ad_mat[, -1, drop = FALSE])
      dp_total   <- ref_counts + alt_counts
      vafs       <- alt_counts / dp_total
      dp_vals    <- dp_total  # derived depth
    }
    
    dp_field <- geno(vcf)$DP
    if (!is.null(dp_field)) {
      dp_vals <- as.numeric(dp_field[, 1])
    } else if (!exists("dp_vals")) {
      dp_vals <- rep(NA_real_, length(vafs))
    }
    
    # ── Depth filter (variant DP >= 30) ──────────────────────────────────────
    pass_depth <- !is.na(dp_vals) & dp_vals >= 30
    if (!any(pass_depth)) return(NULL)
    vafs    <- vafs[pass_depth]
    dp_vals <- dp_vals[pass_depth]
    pos_vec <- start(rowRanges(vcf))[pass_depth]
    # ─────────────────────────────────────────────────────────────────────────
    
    tibble(
      vaf = as.numeric(vafs),
      dp  = as.numeric(dp_vals),
      pos = pos_vec
    )
    
  }, error = function(e) {
    warning(sprintf("VCF read error at %s:%d-%d — %s", chrom, start, end, e$message))
    NULL
  })
}

# ============================================================
# EMPTY TIBBLE TEMPLATE — safe fallback when no variants found
# ============================================================

empty_vaf_tibble <- function() {
  tibble(
    vaf               = numeric(),
    dp                = numeric(),
    pos               = integer(),
    Sample            = character(),
    sample_id_trimmed = character(),
    run               = character(),
    Gene              = character(),
    CNV_Type          = character(),
    Chromosome        = character(),
    CNV_Start         = integer(),
    CNV_End           = integer(),
    acgh_status       = character(),
    vaf_status        = character()
  )
}

# ============================================================
# GENERIC LOOP HELPER: EXTRACT VAFs FOR A SET OF CNV CALLS
# ============================================================

extract_vafs_for_calls <- function(calls_df, vcf_lookup, vaf_upper = 0.95) {
  results <- lapply(seq_len(nrow(calls_df)), function(i) {
    
    row      <- calls_df[i, ]
    dns_id   <- row$sample_id_trimmed
    vcf_path <- vcf_lookup[dns_id]
    
    if (is.na(vcf_path) || !file.exists(vcf_path)) return(NULL)
    
    vafs_df <- extract_vafs_in_region(
      vcf_path = vcf_path,
      chrom    = row$Chromosome,
      start    = row$Start,
      end      = row$End
    )
    
    if (is.null(vafs_df) || nrow(vafs_df) == 0) return(NULL)
    
    vafs_df %>%
      filter(!is.na(vaf), !is.na(dp)) %>%
      filter(vaf < vaf_upper) %>%
      mutate(
        Sample            = row$Sample,
        sample_id_trimmed = dns_id,
        run               = row$run,
        Gene              = row$Gene,
        CNV_Type          = row$CNV_Type,
        Chromosome        = row$Chromosome,
        CNV_Start         = row$Start,
        CNV_End           = row$End,
        acgh_status       = row$acgh_status
      )
  })
  
  out <- bind_rows(results)
  if (ncol(out) == 0) return(empty_vaf_tibble())
  out
}

# ============================================================
# ANNOTATE VAF STATUS
# Adds vaf_status column: "concordant" or "discordant"
# Duplication: concordant if CI overlaps 1/3 or 2/3
# Deletion:    concordant if CI does NOT overlap 0.5
#              (discordant = retained heterozygosity)
# ============================================================

annotate_vaf_status_dup <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      vaf_status = case_when(
        # CI includes 0.5 AND (1/3 or 2/3) → ambiguous, neutral
        vaf_within_ci_het(vaf, dp, target = 0.5) &
          vaf_within_ci(vaf, dp, targets = c(1/3, 2/3))     ~ "neutral",
        # CI includes 1/3 or 2/3 but NOT 0.5 → supportive of duplication
        vaf_within_ci(vaf, dp, targets = c(1/3, 2/3))       ~ "concordant",
        # CI excludes 1/3, 2/3, and 0.5 → not supportive
        TRUE                                                 ~ "discordant"
      )
    ) %>%
    ungroup()
}


annotate_vaf_status_del <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      vaf_status = if_else(
        !vaf_within_ci_het(vaf, dp, target = 0.0) & 
          !vaf_within_ci_het(vaf, dp, target = 1.0),
        "discordant",   # CI excludes both 0 and 1.0 = indicative of heterozygosity = discordant with deletion-associated homozygosity
        "neutral"       # CI includes 0 or 1.0 = consistent with LOH/homozygous = concordant with deletion
      )
    ) %>%
    ungroup()
}

# ============================================================
# IDENTIFY CALL SETS
# ============================================================

discordant_dup_calls <- decon_with_acgh %>%
  filter(acgh_status == "discordant", CNV_Type == "duplication")

discordant_del_calls <- decon_with_acgh %>%
  filter(acgh_status == "discordant", CNV_Type == "deletion")

concordant_dup_calls <- decon_with_acgh %>%
  filter(acgh_status == "concordant", CNV_Type == "duplication")

concordant_del_calls <- decon_with_acgh %>%
  filter(acgh_status == "concordant", CNV_Type == "deletion")

# ============================================================
# EXTRACT AND ANNOTATE VAFs — DUPLICATIONS (FP + TP)
# ============================================================

vaf_all_dup_fp <- extract_vafs_for_calls(discordant_dup_calls, vcf_lookup, vaf_upper = 0.95) %>%
  annotate_vaf_status_dup()

vaf_all_dup_tp <- extract_vafs_for_calls(concordant_dup_calls, vcf_lookup, vaf_upper = 0.95) %>%
  annotate_vaf_status_dup()

# ============================================================
# EXTRACT AND ANNOTATE VAFs — DELETIONS (FP + TP)
# ============================================================

vaf_all_del_fp <- extract_vafs_for_calls(discordant_del_calls, vcf_lookup, vaf_upper = 0.95) %>%
  filter(vaf > 0) %>%
  annotate_vaf_status_del()

vaf_all_del_tp <- extract_vafs_for_calls(concordant_del_calls, vcf_lookup, vaf_upper = 0.95) %>%
  filter(vaf > 0) %>%
  annotate_vaf_status_del()

# ============================================================
# SUMMARISE — FP DUPLICATIONS
# ============================================================

fp_dup_calls_with_vaf <- vaf_all_dup_fp %>%
  distinct(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type)

fp_vaf_any_discordant_dup <- vaf_all_dup_fp %>%
  group_by(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type) %>%
  summarise(
    n_vafs          = n(),
    n_discordant    = sum(vaf_status == "discordant", na.rm = TRUE),
    prop_discordant = n_discordant / n_vafs,
    .groups = "drop"
  ) %>%
  filter(n_discordant >= 1)

fp_vaf_majority_discordant_dup <- fp_vaf_any_discordant_dup %>%
  filter(prop_discordant >= 0.5)

# ============================================================
# SUMMARISE — TP DUPLICATIONS
# ============================================================

tp_dup_calls_with_vaf <- vaf_all_dup_tp %>%
  distinct(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type)

tp_vaf_any_discordant_dup <- vaf_all_dup_tp %>%
  group_by(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type) %>%
  summarise(
    n_vafs          = n(),
    n_discordant    = sum(vaf_status == "discordant", na.rm = TRUE),
    prop_discordant = n_discordant / n_vafs,
    .groups = "drop"
  ) %>%
  filter(n_discordant >= 1)

tp_vaf_majority_discordant_dup <- tp_vaf_any_discordant_dup %>%
  filter(prop_discordant >= 0.5)

# ============================================================
# SUMMARISE — FP DELETIONS
# ============================================================

fp_del_calls_with_vaf <- vaf_all_del_fp %>%
  distinct(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type)

fp_vaf_discordant_del <- vaf_all_del_fp %>%
  group_by(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type) %>%
  summarise(
    n_vafs             = n(),
    any_discordant_vaf = any(vaf_status == "discordant", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(any_discordant_vaf)

fp_vaf_any_discordant_del <- vaf_all_del_fp %>%
  group_by(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type) %>%
  summarise(
    n_vafs          = n(),
    n_discordant    = sum(vaf_status == "discordant", na.rm = TRUE),
    prop_discordant = n_discordant / n_vafs,
    .groups = "drop"
  ) %>%
  filter(n_discordant >= 1)

fp_vaf_majority_discordant_del <- fp_vaf_any_discordant_del %>%
  filter(prop_discordant >= 0.5)

# ============================================================
# SUMMARISE — TP DELETIONS
# ============================================================

tp_del_calls_with_vaf <- vaf_all_del_tp %>%
  distinct(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type)

tp_vaf_discordant_del <- vaf_all_del_tp %>%
  group_by(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type) %>%
  summarise(
    n_vafs             = n(),
    any_discordant_vaf = any(vaf_status == "discordant", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(any_discordant_vaf)

tp_vaf_any_discordant_del <- vaf_all_del_tp %>%
  group_by(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type) %>%
  summarise(
    n_vafs          = n(),
    n_discordant    = sum(vaf_status == "discordant", na.rm = TRUE),
    prop_discordant = n_discordant / n_vafs,
    .groups = "drop"
  ) %>%
  filter(n_discordant >= 1)

tp_vaf_majority_discordant_del <- tp_vaf_any_discordant_del %>%
  filter(prop_discordant >= 0.5)

# ============================================================
# PRINT SUMMARY
# ============================================================

cat("=== DUPLICATION CALLS ===\n")
cat(sprintf("FP duplication calls with ≥1 informative SNV VAF (< 0.95):      %d / %d\n",
            nrow(fp_dup_calls_with_vaf), nrow(discordant_dup_calls)))
cat(sprintf("FP duplication calls with ≥1 VAF discordant from 1/3 or 2/3:   %d / %d\n",
            nrow(fp_vaf_any_discordant_dup), nrow(fp_dup_calls_with_vaf)))
cat(sprintf("FP duplication calls with ≥50%% discordant VAFs:                 %d / %d\n",
            nrow(fp_vaf_majority_discordant_dup), nrow(fp_dup_calls_with_vaf)))
cat(sprintf("TP duplication calls with ≥1 informative SNV VAF (< 0.95):      %d / %d\n",
            nrow(tp_dup_calls_with_vaf), nrow(concordant_dup_calls)))
cat(sprintf("TP duplication calls with ≥1 VAF discordant from 1/3 or 2/3:   %d / %d\n",
            nrow(tp_vaf_any_discordant_dup), nrow(tp_dup_calls_with_vaf)))
cat(sprintf("TP duplication calls with ≥50%% discordant VAFs:                 %d / %d\n",
            nrow(tp_vaf_majority_discordant_dup), nrow(tp_dup_calls_with_vaf)))

cat("\n=== DELETION CALLS ===\n")
cat(sprintf("FP deletion calls with ≥1 informative SNV VAF (> 0 and < 1):    %d / %d\n",
            nrow(fp_del_calls_with_vaf), nrow(discordant_del_calls)))
cat(sprintf("FP deletion calls with ≥1 discordant VAF:              %d / %d\n",
            nrow(fp_vaf_discordant_del), nrow(fp_del_calls_with_vaf)))
cat(sprintf("FP deletion calls with ≥50%% discordant VAFs:                  %d / %d\n",
            nrow(fp_vaf_majority_discordant_del), nrow(fp_del_calls_with_vaf)))
cat(sprintf("TP deletion calls with ≥1 informative SNV VAF (> 0 and < 1):    %d / %d\n",
            nrow(tp_del_calls_with_vaf), nrow(concordant_del_calls)))
cat(sprintf("TP deletion calls with ≥1 VAF consistent with 0.5 (unexpected): %d / %d\n",
            nrow(tp_vaf_discordant_del), nrow(tp_del_calls_with_vaf)))
cat(sprintf("TP deletion calls with ≥50%% discordant VAFs (unexpected):     %d / %d\n",
            nrow(tp_vaf_majority_discordant_del), nrow(tp_del_calls_with_vaf)))


# ============================================================
# OPTIMISED FRAMEWORK NUMBERS
# ============================================================
# ============================================================
# STAGE 1: REMOVE OUTLIER SAMPLES
# ============================================================

outlier_samples <- all_samples %>%
  filter(outlier_sample == "YES") %>%
  dplyr::select(sample_id_trimmed, run)

filter_outlier_samples <- function(df, outlier_df) {
  df %>% anti_join(outlier_df, by = c("sample_id_trimmed", "run"))
}

decon_no_outlier    <- filter_outlier_samples(decon_with_acgh, outlier_samples)
clearcnv_no_outlier <- filter_outlier_samples(clearcnv_with_acgh, outlier_samples)

stage1_table <- dplyr::bind_rows(
  decon_no_outlier    %>% dplyr::mutate(caller = "DECoN"),
  clearcnv_no_outlier %>% dplyr::mutate(caller = "clearCNV")
) %>%
  dplyr::mutate(acgh_category = case_when(
    acgh_status %in% c("tbc", "inconclusive") ~ "tbc/inconclusive",
    TRUE ~ acgh_status
  )) %>%
  dplyr::count(caller, acgh_category) %>%
  arrange(caller, acgh_category)

print("=== Stage 1: Calls after outlier removal ===")
print(stage1_table)

# ============================================================
# HELPER: MEAN READ COUNT PER CNV CALL
#         - Uses exons overlapping CNV coordinates
#         - Works for all aCGH statuses
# ============================================================

get_cnv_region_read_count <- function(calls_df) {
  
  runs_needed <- unique(calls_df$run)
  result_list <- list()
  
  for (run_num in runs_needed) {
    rdata_path <- sprintf(
      "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/exome_depth/RData/run%d_bams_9genes_25bp.fix.sorted.RData",
      run_num
    )
    
    if (!file.exists(rdata_path)) {
      warning(sprintf("Run %d: RData not found, skipping.", run_num))
      next
    }
    
    load(rdata_path)  # loads my.count
    
    my.count.dafr <- as(my.count, "data.frame")
    sample_cols   <- colnames(my.count.dafr)[6:ncol(my.count.dafr)]
    count_matrix  <- my.count.dafr[, sample_cols] %>% mutate_all(as.numeric)
    
    colnames(count_matrix) <- colnames(count_matrix) %>%
      str_remove("\\.bam$") %>%
      str_remove("\\.hq\\.sorted\\.marked.*$") %>%
      str_replace_all("\\.", "-")
    
    exon_coords <- data.frame(
      chromosome = as.character(gsub("^chr", "", my.count.dafr[["chromosome"]])),
      start      = my.count.dafr[["start"]],
      end        = my.count.dafr[["end"]],
      exon_idx   = seq_len(nrow(my.count.dafr))
    )
    
    run_calls <- calls_df %>%
      filter(run == run_num) %>%
      mutate(
        Sample_normalised = str_remove(Sample, "\\.hq\\.sorted\\.marked.*$")
      )
    
    row_results <- lapply(seq_len(nrow(run_calls)), function(i) {
      cnv        <- run_calls[i, ]
      cnv_chr    <- as.character(cnv$Chromosome)
      cnv_start  <- cnv$Start
      cnv_end    <- cnv$End
      cnv_sample <- cnv$Sample_normalised
      
      overlapping_exons <- exon_coords %>%
        filter(
          chromosome == cnv_chr,
          start <= cnv_end,
          end >= cnv_start
        ) %>%
        pull(exon_idx)
      
      if (length(overlapping_exons) == 0 || !cnv_sample %in% colnames(count_matrix)) {
        mean_rc <- NA_real_
      } else {
        exon_counts <- count_matrix[overlapping_exons, cnv_sample]
        mean_rc <- mean(exon_counts, na.rm = TRUE)
      }
      
      cnv %>%
        mutate(mean_read_count = mean_rc) %>%
        dplyr::select(-Sample_normalised)
    })
    
    result_list[[as.character(run_num)]] <- bind_rows(row_results)
  }
  
  bind_rows(result_list)
}

# ============================================================
# ADD MEAN READ COUNT TO STAGE 1 CALLS
# ============================================================

decon_stage1_rc    <- get_cnv_region_read_count(decon_no_outlier)
clearcnv_stage1_rc <- get_cnv_region_read_count(clearcnv_no_outlier)

# ============================================================
# STAGE 2: APPLY FILTERS
#          DECoN: BF >= 10 and mean_read_count >= 100
#          clearCNV: mean_read_count >= 100 only
# ============================================================

decon_filtered <- decon_stage1_rc %>%
  filter(Quality >= 10, mean_read_count >= 100)

clearcnv_filtered <- clearcnv_stage1_rc %>%
  filter(mean_read_count >= 100)

stage2_table <- bind_rows(
  decon_filtered    %>% mutate(caller = "DECoN"),
  clearcnv_filtered %>% mutate(caller = "clearCNV")
) %>%
  mutate(acgh_category = case_when(
    acgh_status %in% c("tbc", "inconclusive") ~ "tbc/inconclusive",
    TRUE ~ acgh_status
  )) %>%
  dplyr::count(caller, acgh_category) %>%
  arrange(caller, acgh_category)

print("=== Stage 2: Calls after filters ===")
print(stage2_table)

# ============================================================
# STAGE 3: PRIORITISATION MATRIX
#          - high_confidence: called by BOTH DECoN AND clearCNV
#          - cnv_need_manual_inspection: called by EITHER caller only
# ============================================================

find_reciprocal_overlaps <- function(decon_df, clearcnv_df,
                                     reciprocal_threshold = 0.5) {
  
  make_gr <- function(df, caller_label) {
    GRanges(
      seqnames = df$Chromosome,
      ranges   = IRanges(start = df$Start, end = df$End),
      Sample   = df$Sample,
      sample_id_trimmed = df$sample_id_trimmed,
      run      = df$run,
      Gene     = df$Gene,
      CNV_Type = df$CNV_Type,
      Quality  = df$Quality,
      acgh_status        = df$acgh_status,
      cyto_overlap       = df$cyto_overlap,
      mean_read_count    = df$mean_read_count,
      caller   = caller_label
    )
  }
  
  gr_decon    <- make_gr(decon_df,    "DECoN")
  gr_clearcnv <- make_gr(clearcnv_df, "clearCNV")
  
  # Find all pairwise overlaps
  hits <- findOverlaps(gr_decon, gr_clearcnv, ignore.strand = TRUE)
  
  # Keep only same-sample, same-CNV_Type, same-run overlaps
  same_sample <- gr_decon$Sample[queryHits(hits)]   == gr_clearcnv$Sample[subjectHits(hits)]
  same_type   <- gr_decon$CNV_Type[queryHits(hits)] == gr_clearcnv$CNV_Type[subjectHits(hits)]
  same_run    <- gr_decon$run[queryHits(hits)]       == gr_clearcnv$run[subjectHits(hits)]
  
  valid_hits <- hits[same_sample & same_type & same_run]
  
  # Compute reciprocal overlap fraction
  decon_ranges    <- ranges(gr_decon)[queryHits(valid_hits)]
  clearcnv_ranges <- ranges(gr_clearcnv)[subjectHits(valid_hits)]
  
  intersection_width <- width(pintersect(decon_ranges, clearcnv_ranges))
  decon_width        <- width(decon_ranges)
  clearcnv_width     <- width(clearcnv_ranges)
  
  reciprocal_overlap <- pmin(
    intersection_width / decon_width,
    intersection_width / clearcnv_width
  )
  
  # Flag pairs meeting reciprocal threshold
  concordant_hits <- valid_hits[reciprocal_overlap >= reciprocal_threshold]
  
  decon_concordant_idx    <- unique(queryHits(concordant_hits))
  clearcnv_concordant_idx <- unique(subjectHits(concordant_hits))
  
  # Build high_confidence set (union of matched DECoN + clearCNV calls)
  # Anchor on DECoN row; attach best clearCNV match metadata
  best_match <- as.data.frame(concordant_hits) %>%
    mutate(recip_overlap = reciprocal_overlap[reciprocal_overlap >= reciprocal_threshold]) %>%
    group_by(queryHits) %>%
    slice_max(recip_overlap, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  high_confidence <- as.data.frame(gr_decon[best_match$queryHits]) %>%
    rename(Start = start, End = end, Chromosome = seqnames) %>%
    dplyr::select(-strand, -width) %>%
    mutate(
      priority           = "high_confidence (called by DECoN and ClearCNV)",
      clearcnv_quality   = gr_clearcnv$Quality[best_match$subjectHits],
      reciprocal_overlap = best_match$recip_overlap
    )
  
  # DECoN-only calls
  decon_only <- as.data.frame(gr_decon[-decon_concordant_idx]) %>%
    rename(Start = start, End = end, Chromosome = seqnames) %>%
    dplyr::select(-strand, -width) %>%
    mutate(
      priority           = "cnv_need_manual_inspection",
      clearcnv_quality   = NA_real_,
      reciprocal_overlap = NA_real_
    )
  
  # clearCNV-only calls
  clearcnv_only <- as.data.frame(gr_clearcnv[-clearcnv_concordant_idx]) %>%
    rename(Start = start, End = end, Chromosome = seqnames) %>%
    dplyr::select(-strand, -width) %>%
    mutate(
      priority           = "cnv_need_manual_inspection",
      clearcnv_quality   = Quality,
      reciprocal_overlap = NA_real_
    )
  
  bind_rows(high_confidence, decon_only, clearcnv_only)
}

# ============================================================
# RUN STAGE 3
# ============================================================

stage3_calls <- find_reciprocal_overlaps(
  decon_df    = decon_filtered,
  clearcnv_df = clearcnv_filtered,
  reciprocal_threshold = 0.5
)

# ============================================================
# STAGE 3 SUMMARY TABLE
# ============================================================

stage3_table <- stage3_calls %>%
  mutate(acgh_category = case_when(
    acgh_status %in% c("tbc", "inconclusive") ~ "tbc/inconclusive",
    TRUE ~ acgh_status
  )) %>%
  dplyr::count(caller, priority, acgh_category) %>%
  arrange(caller, priority, acgh_category)

print("=== Stage 3: Prioritisation matrix ===")
print(stage3_table)

manual_inspection_calls <- stage3_calls %>%
  filter(priority == "cnv_need_manual_inspection") %>%
  mutate(Chromosome = as.numeric(gsub("^chr", "", as.character(Chromosome))))

##### Check VAF concordance for manually inspected calls! ########
# ============================================================
# APPLY SAME ANALYSIS TO MANUAL INSPECTION CALLS
# ============================================================

manual_dup_calls <- manual_inspection_calls %>%
  filter(CNV_Type == "duplication")

manual_del_calls <- manual_inspection_calls %>%
  filter(CNV_Type == "deletion")

# Extract and annotate duplication VAFs
manual_vaf_all_dup <- extract_vafs_for_calls(manual_dup_calls, vcf_lookup, vaf_upper = 0.95) %>%
  annotate_vaf_status_dup()

# Extract and annotate deletion VAFs
manual_vaf_all_del <- extract_vafs_for_calls(manual_del_calls, vcf_lookup, vaf_upper = 0.95) %>%
  filter(vaf > 0) %>%
  annotate_vaf_status_del()

# Summarise duplication calls
manual_dup_calls_with_vaf <- manual_vaf_all_dup %>%
  distinct(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type)

manual_dup_summary <- manual_vaf_all_dup %>%
  group_by(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type, acgh_status) %>%
  summarise(
    n_vafs          = n(),
    n_concordant    = sum(vaf_status == "concordant", na.rm = TRUE),
    n_discordant    = sum(vaf_status == "discordant", na.rm = TRUE),
    prop_discordant = n_discordant / n_vafs,
    any_discordant  = any(vaf_status == "discordant", na.rm = TRUE),
    .groups = "drop"
  )

# Summarise deletion calls
manual_del_calls_with_vaf <- manual_vaf_all_del %>%
  distinct(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type)

manual_del_summary <- manual_vaf_all_del %>%
  group_by(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type, acgh_status) %>%
  summarise(
    n_vafs          = n(),
    n_concordant    = sum(vaf_status == "concordant", na.rm = TRUE),
    n_discordant    = sum(vaf_status == "discordant", na.rm = TRUE),
    prop_discordant = n_discordant / n_vafs,
    any_discordant  = any(vaf_status == "discordant", na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# STAGE 4: HIGH-CONFIDENCE CLASSIFICATION
# ============================================================

# ── 4a. Split stage3 calls by priority ───────────────────────────────────────

# Normalise Chromosome to character for both branches before splitting
# (avoids factor vs double mismatch at bind_rows later)
stage3_calls <- stage3_calls %>%
  mutate(Chromosome = as.character(gsub("^chr", "", as.character(Chromosome))))

# Calls agreed upon by both callers → automatically high-confidence
dual_caller_hc <- stage3_calls %>%
  filter(priority == "high_confidence (called by DECoN and ClearCNV)") %>%
  mutate(hc_classification = "high_confidence")

# Calls from one caller only → require RR + VAF manual inspection criteria
manual_inspection_calls <- stage3_calls %>%
  filter(priority == "cnv_need_manual_inspection")

# ── 4b. Load BED file (6-col, 0-based coords, chr-prefixed) ──────────────────

bed_file <- read.table(
  "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/exome_depth/bed_files/9genes_25bp.fix.sorted.bed",
  header = FALSE, sep = "\t",
  col.names = c("chr", "start", "end", "name", "score", "strand")
) %>%
  mutate(
    gene  = str_extract(name, "^[^_]+"),
    start = start + 1L
  )

bed_gr <- GRanges(
  seqnames = bed_file$chr,
  ranges   = IRanges(start = bed_file$start, end = bed_file$end),
  gene     = bed_file$gene
)

# ── 4c. Count exons per CNV call (manual inspection calls only) ───────────────

count_exons_in_call <- function(calls_df) {
  calls_gr <- GRanges(
    seqnames = paste0("chr", gsub("^chr", "", as.character(calls_df$Chromosome))),
    ranges   = IRanges(start = calls_df$Start, end = calls_df$End)
  )
  
  hits <- findOverlaps(calls_gr, bed_gr, ignore.strand = TRUE)
  
  hit_df <- data.frame(
    call_idx = queryHits(hits),
    bed_gene = bed_gr$gene[subjectHits(hits)]
  ) %>%
    mutate(cnv_gene = calls_df$Gene[call_idx]) %>%
    filter(bed_gene == cnv_gene) %>%
    dplyr::count(call_idx, name = "n_exons_in_call")
  
  calls_df %>%
    mutate(call_idx = row_number()) %>%
    left_join(hit_df, by = "call_idx") %>%
    mutate(n_exons_in_call = replace_na(n_exons_in_call, 0L)) %>%
    dplyr::select(-call_idx)
}

manual_with_exons <- count_exons_in_call(manual_inspection_calls)

# ── 4d. Compute per-call RR features (manual inspection calls only) ───────────

rr_input <- manual_with_exons %>%
  mutate(
    acgh_status = if_else(
      acgh_status %in% c("tbc", "inconclusive"),
      "unknown",
      replace_na(acgh_status, "unknown")
    )
  )

all_rr_exons <- get_decon_exon_read_ratios(rr_input)

rr_pass <- all_rr_exons %>%
  filter(exon_type == "in_call") %>%
  mutate(expected_rr = if_else(CNV_Type == "deletion", 0.5, 1.5)) %>%
  group_by(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type) %>%
  summarise(
    n_in_call_exons = n(),
    rr1_pass  = any(abs(read_ratio - expected_rr) <= 0.15, na.rm = TRUE),
    rr2_pass  = if_else(
      CNV_Type[1L] == "deletion",
      !any(read_ratio >= 1.35, na.rm = TRUE),
      !any(read_ratio <= 0.65, na.rm = TRUE)
    ),
    rr_spread = max(read_ratio, na.rm = TRUE) - min(read_ratio, na.rm = TRUE),
    rr3_pass  = if_else(n_in_call_exons > 1L, rr_spread <= 0.5, TRUE),
    rr_pass_all = rr1_pass & rr2_pass & rr3_pass,
    .groups = "drop"
  )

# ── 4e. Extract and annotate VAFs (manual inspection calls only) ──────────────

all_vaf_del <- extract_vafs_for_calls(
  manual_with_exons %>% filter(CNV_Type == "deletion"),
  vcf_lookup, vaf_upper = 0.95
) %>%
  filter(vaf > 0) %>%
  annotate_vaf_status_del()

del_vaf_summary <- all_vaf_del %>%
  group_by(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type) %>%
  summarise(
    n_informative_vafs = n(),
    n_discordant       = sum(vaf_status == "discordant", na.rm = TRUE),
    any_discordant     = any(vaf_status == "discordant", na.rm = TRUE),
    .groups = "drop"
  )

all_vaf_dup <- extract_vafs_for_calls(
  manual_with_exons %>% filter(CNV_Type == "duplication"),
  vcf_lookup, vaf_upper = 0.95
) %>%
  annotate_vaf_status_dup()

dup_vaf_summary <- all_vaf_dup %>%
  group_by(Sample, Gene, Chromosome, CNV_Start, CNV_End, CNV_Type) %>%
  summarise(
    n_informative_vafs = n(),
    n_concordant       = sum(vaf_status == "concordant", na.rm = TRUE),
    n_discordant       = sum(vaf_status == "discordant", na.rm = TRUE),
    prop_concordant    = n_concordant / n_informative_vafs,
    any_discordant     = any(vaf_status == "discordant", na.rm = TRUE),
    .groups = "drop"
  )

# ── 4f. Join RR + VAF features onto manual inspection calls ──────────────────

stage4_manual <- manual_with_exons %>%
  left_join(
    rr_pass %>% rename(Start = CNV_Start, End = CNV_End),
    by = c("Sample", "Gene", "Chromosome", "Start", "End", "CNV_Type")
  ) %>%
  left_join(
    del_vaf_summary %>%
      rename(Start = CNV_Start, End = CNV_End) %>%
      rename(del_n_informative_vafs = n_informative_vafs,
             del_n_discordant       = n_discordant,
             del_any_discordant     = any_discordant),
    by = c("Sample", "Gene", "Chromosome", "Start", "End", "CNV_Type")
  ) %>%
  left_join(
    dup_vaf_summary %>%
      rename(Start = CNV_Start, End = CNV_End) %>%
      rename(dup_n_informative_vafs = n_informative_vafs,
             dup_n_concordant       = n_concordant,
             dup_n_discordant       = n_discordant,
             dup_prop_concordant    = prop_concordant,
             dup_any_discordant     = any_discordant),
    by = c("Sample", "Gene", "Chromosome", "Start", "End", "CNV_Type")
  ) %>%
  mutate(
    n_informative_vafs = if_else(CNV_Type == "deletion",
                                 replace_na(del_n_informative_vafs, 0L),
                                 replace_na(dup_n_informative_vafs, 0L)),
    n_discordant       = if_else(CNV_Type == "deletion",
                                 replace_na(del_n_discordant, 0L),
                                 replace_na(dup_n_discordant, 0L)),
    any_discordant     = if_else(CNV_Type == "deletion",
                                 replace_na(del_any_discordant, FALSE),
                                 replace_na(dup_any_discordant, FALSE)),
    n_concordant       = replace_na(dup_n_concordant,    0L),
    prop_concordant    = replace_na(dup_prop_concordant, 0)
  ) %>%
  dplyr::select(-starts_with("del_"), -starts_with("dup_"))

# ── 4g. Apply high-confidence decision rules to manual inspection calls ───────

stage4_manual <- stage4_manual %>%
  mutate(
    is_single_exon = replace_na(n_exons_in_call == 1L, FALSE),
    
    vaf_hc = case_when(
      
      # ── DELETION ─────────────────────────────────────────────────────────────
      CNV_Type == "deletion" & any_discordant                              ~ FALSE,
      CNV_Type == "deletion" & !is_single_exon                            ~ FALSE,
      CNV_Type == "deletion" & is_single_exon & n_informative_vafs == 0L  ~ TRUE,
      CNV_Type == "deletion" & is_single_exon & n_informative_vafs > 0 & !any_discordant  ~ TRUE,
      
      # ── DUPLICATION ──────────────────────────────────────────────────────────
      CNV_Type == "duplication" & !is_single_exon & n_informative_vafs == 0L          ~ FALSE,
      CNV_Type == "duplication" & !is_single_exon & n_informative_vafs > 0 & prop_concordant > 0.5 ~ TRUE,
      CNV_Type == "duplication" & is_single_exon & n_informative_vafs == 0L           ~ TRUE,
      CNV_Type == "duplication" & is_single_exon & n_informative_vafs > 0 & prop_concordant > 0.5 & n_discordant <= 1L ~ TRUE
    ),
    
    hc_classification = case_when(
      !replace_na(rr_pass_all, FALSE)  ~ "low_confidence_rr_fail",
      isFALSE(vaf_hc)                  ~ "low_confidence_vaf_fail",
      vaf_hc == TRUE                   ~ "high_confidence",   # rr_pass_all TRUE is already guaranteed by row 1
      is.na(vaf_hc)                    ~ "uncertain_no_vaf",
      TRUE                             ~ "uncertain"
    )
  )

# ── 4h. Combine dual-caller HC calls with manual inspection results ───────────

stage4_calls <- bind_rows(
  dual_caller_hc,
  stage4_manual
)

# ── 4i. Isolate all high-confidence calls ────────────────────────────────────

high_confidence_calls <- stage4_calls %>%
  filter(hc_classification == "high_confidence")

# ── 4j. Summary ──────────────────────────────────────────────────────────────

stage4_summary <- stage4_calls %>%
  mutate(acgh_category = case_when(
    acgh_status %in% c("tbc", "inconclusive") ~ "tbc/inconclusive",
    TRUE ~ acgh_status
  )) %>%
  dplyr::count(CNV_Type, hc_classification, acgh_category) %>%
  arrange(CNV_Type, hc_classification, acgh_category)

print("=== Stage 4: High-confidence classification ===")
print(stage4_summary)

cat(sprintf("\nTotal high-confidence calls: %d\n", nrow(high_confidence_calls)))
cat(sprintf("  From dual-caller agreement: %d\n", nrow(dual_caller_hc)))
cat(sprintf("  From manual inspection:     %d\n",
            sum(stage4_manual$hc_classification == "high_confidence", na.rm = TRUE)))
cat(sprintf("  Deletions:    %d\n", sum(high_confidence_calls$CNV_Type == "deletion")))
cat(sprintf("  Duplications: %d\n", sum(high_confidence_calls$CNV_Type == "duplication")))

# write.csv(
#   high_confidence_calls,
#   file = "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_2/Optimized_CNV_Detection_Pipeline/rr_of_12_highconf_calls_with_vaf/high_confidence_calls.csv",
#   row.names = FALSE
# )

# ============================================================
# SUMMARY TABLE: MEAN READ COUNT FOR CONCORDANT/DISCORDANT
# FROM STAGE 1
# ============================================================

decon_rc <- decon_stage1_rc %>%
  filter(acgh_status %in% c("concordant", "discordant"))

clearcnv_rc <- clearcnv_stage1_rc %>%
  filter(acgh_status %in% c("concordant", "discordant"))

rc_summary_table <- bind_rows(
  decon_rc    %>% mutate(caller = "DECoN"),
  clearcnv_rc %>% mutate(caller = "clearCNV")
) %>%
  group_by(caller, acgh_status) %>%
  summarise(
    n_calls              = n(),
    mean_read_count      = round(mean(mean_read_count, na.rm = TRUE), 2),
    median_read_count    = round(median(mean_read_count, na.rm = TRUE), 2),
    n_missing_read_count = sum(is.na(mean_read_count)),
    .groups = "drop"
  ) %>%
  arrange(caller, acgh_status)

print("=== Mean read count for concordant/discordant calls (Stage 1) ===")
print(rc_summary_table)


# ============================================================
# Plotting the rc of concordant vs discordant calls
# ============================================================
decon_rc_plot <- decon_rc %>%
  filter(acgh_status %in% c("concordant", "discordant")) %>%
  mutate(
    call_status = recode(acgh_status,
                         "concordant" = "True Positive",
                         "discordant" = "False Positive"),
    call_status = factor(call_status, levels = c("True Positive", "False Positive"))
  )

status_colours <- c("True Positive" = "#2196F3", "False Positive" = "#F44336")

p_violin <- ggplot(decon_rc_plot, aes(x = call_status, y = mean_read_count, fill = call_status)) +
  geom_violin(alpha = 0.5, trim = FALSE, colour = "grey40") +
  geom_jitter(aes(colour = call_status), width = 0.15, size = 1.8, alpha = 0.8) +
  scale_fill_manual(values   = status_colours) +
  scale_colour_manual(values = status_colours) +
  scale_y_continuous(labels  = scales::comma) +
  labs(
    title = "Mean Exon-level Read Depth of True- and False-Positive Per-Run CNV Calls from DECoN and ClearCNV",
    x     = "Call Classification",
    y     = "Mean Read Depth (Reads)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_violin)

# ggsave(
#   "cnv_read_count_TP_FP_violin_jitter.jpeg",
#   plot   = p_violin,
#   width  = 6,
#   height = 6,
#   units  = "in",
#   dpi    = 300
# )



# ============================================================
# OPTIONAL: WRITE ALL TABLES
# ============================================================

# write.csv(stage1_table,      "stage1_post_outlier_removal.csv",          row.names = FALSE)
# write.csv(stage2_table,      "stage2_post_quality_filter.csv",           row.names = FALSE)
# write.csv(rc_summary_table,  "stage3_read_count_concordant_discordant.csv", row.names = FALSE)

# Full per-call table with read counts (useful for inspection)
# write.csv(
#   bind_rows(decon_rc %>% mutate(caller = "DECoN"),
#             clearcnv_rc %>% mutate(caller = "clearCNV")),
#   "cnv_calls_with_region_read_counts.csv",
#   row.names = FALSE
# )
