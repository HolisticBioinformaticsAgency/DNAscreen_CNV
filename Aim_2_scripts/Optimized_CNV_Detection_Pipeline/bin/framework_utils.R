# ============================================================
# Shared functions for the optimised CNV framework
# (ported from optimized_framework_calls.R / plot_rr_with_vaf.R)
# Sourced by the other bin/*.R scripts.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(GenomicRanges)
  library(IRanges)
  library(Rsamtools)
})

# ============================================================
# NAMES / IDS
# ============================================================

# BAM path / DECoN sample / clearCNV sample -> "DNS-XTHS-0000-A00-DNS000000_S1"
std_sample <- function(x) {
  x %>%
    basename() %>%
    str_remove("\\.bam$") %>%
    str_remove("\\.hq\\.sorted\\.marked.*$")
}

std_chr <- function(x) gsub("^chr", "", as.character(x))

make_call_id <- function(run, caller, sample, chr, start, end, cnv_type) {
  sprintf("%s_%s_%s_chr%s:%s-%s_%s", run, caller, sample, std_chr(chr),
          format(start, scientific = FALSE, trim = TRUE),
          format(end,   scientific = FALSE, trim = TRUE), cnv_type)
}

# Self-contained HTML needs pandoc (module load pandoc); otherwise writes <name>_files/ alongside
save_widget_html <- function(widget, path) {
  htmlwidgets::saveWidget(widget, file = normalizePath(path, mustWork = FALSE),
                          selfcontained = nzchar(Sys.which("pandoc")))
}

write_tsv_safe <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write_tsv(df, path, na = "NA")
}

read_tsv_safe <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) return(tibble())
  suppressMessages(read_tsv(path, col_types = cols(.default = col_guess()),
                            guess_max = 1e5, na = c("", "NA")))
}

# ============================================================
# BED / GENE ANNOTATION
# ============================================================

read_framework_bed <- function(bed_path) {
  read.table(bed_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1:4] %>%
    setNames(c("chr", "start", "end", "name")) %>%
    mutate(
      chr   = std_chr(chr),
      gene  = str_extract(name, "^[^_]+"),
      start = start + 1L                     # BED 0-based -> 1-based
    )
}

gene_strand <- c(BRCA1 = "-", BRCA2 = "+", PALB2 = "-", MSH2 = "+", MSH6 = "+",
                 MLH1 = "+", APOB = "-", PCSK9 = "+", LDLR = "+")

clinically_relevant_transcripts <- c(
  BRCA1 = "NM_007294.4", BRCA2 = "NM_000059.4", PALB2 = "NM_024675.4",
  MSH2  = "NM_000251.3", MSH6  = "NM_000179.3", MLH1  = "NM_000249.4",
  APOB  = "NM_000384.3", PCSK9 = "NM_174936.4", LDLR  = "NM_000527.5"
)

# Number of framework-BED exons of the call's own gene overlapped by the call
count_exons_in_call <- function(calls_df, bed_df) {
  if (nrow(calls_df) == 0) return(calls_df %>% mutate(n_exons_in_call = integer()))
  calls_gr <- GRanges(std_chr(calls_df$Chromosome),
                      IRanges(calls_df$Start, calls_df$End))
  bed_gr   <- GRanges(bed_df$chr, IRanges(bed_df$start, bed_df$end), gene = bed_df$gene)
  hits <- findOverlaps(calls_gr, bed_gr, ignore.strand = TRUE)
  hit_df <- tibble(call_idx = queryHits(hits),
                   bed_gene = bed_gr$gene[subjectHits(hits)]) %>%
    mutate(cnv_gene = calls_df$Gene[call_idx]) %>%
    filter(bed_gene == cnv_gene) %>%
    dplyr::count(call_idx, name = "n_exons_in_call")
  calls_df %>%
    mutate(call_idx = row_number()) %>%
    left_join(hit_df, by = "call_idx") %>%
    mutate(n_exons_in_call = replace_na(n_exons_in_call, 0L)) %>%
    dplyr::select(-call_idx)
}

# ============================================================
# READ COUNTS (ExomeDepth getBamCounts output, as saved by DECoN ReadInBams.R)
# ============================================================

# Returns list(anno = chromosome/start/end/exon/Gene, raw = count matrix,
#              bams = named vector sample -> bam path)
load_counts <- function(rdata_path) {
  env <- new.env()
  load(rdata_path, envir = env)
  df   <- as.data.frame(env$counts)
  samp <- std_sample(env$bams)
  raw  <- as.matrix(df[, 6:ncol(df)])
  storage.mode(raw) <- "numeric"
  raw[is.na(raw)] <- 0
  colnames(raw) <- samp
  anno <- df[, 1:5] %>%
    mutate(chromosome = std_chr(chromosome),
           Gene       = str_extract(as.character(exon), "^[^_]+"))
  list(anno = anno, raw = raw, bams = setNames(env$bams, samp))
}

# Counts-per-thousand normalisation (no GC loess), as in the framework
cpt_normalise <- function(raw) {
  sweep(raw, 2, colSums(raw, na.rm = TRUE), FUN = "/") * 1e3
}

# Average percent difference of each sample vs the median of the other samples
# (used to pick the top-N reference samples)
compute_apd <- function(norm) {
  s_ids <- colnames(norm)
  vapply(s_ids, function(s) {
    ref_median <- apply(norm[, setdiff(s_ids, s), drop = FALSE], 1, median, na.rm = TRUE)
    mean(abs(norm[, s] - ref_median) / (ref_median + 1e-9) * 100, na.rm = TRUE)
  }, numeric(1))
}

top_references <- function(apd, target, n_ref) {
  cand <- setdiff(names(apd), target)
  names(sort(apd[cand]))[seq_len(min(n_ref, length(cand)))]
}

# Mean raw read count over all exons overlapping the call (any gene)
call_mean_read_count <- function(calls_df, counts) {
  if (nrow(calls_df) == 0) return(numeric())
  a <- counts$anno
  vapply(seq_len(nrow(calls_df)), function(i) {
    cnv <- calls_df[i, ]
    idx <- which(a$chromosome == std_chr(cnv$Chromosome) &
                   a$start <= cnv$End & a$end >= cnv$Start)
    if (length(idx) == 0 || !cnv$sample_std %in% colnames(counts$raw)) return(NA_real_)
    mean(counts$raw[idx, cnv$sample_std], na.rm = TRUE)
  }, numeric(1))
}

# Per-exon read ratios for every exon of the call's gene
# (port of get_decon_exon_read_ratios, restricted to retained samples)
exon_read_ratios <- function(calls_df, counts, norm, apd, n_ref = 20) {
  a <- counts$anno
  out <- lapply(seq_len(nrow(calls_df)), function(i) {
    cnv    <- calls_df[i, ]
    target <- cnv$sample_std
    if (!target %in% colnames(norm)) {
      warning(sprintf("Sample %s not found in counts", target))
      return(NULL)
    }
    refs <- top_references(apd, target, n_ref)

    gene_idx <- which(a$Gene == cnv$Gene & a$chromosome == std_chr(cnv$Chromosome))
    if (length(gene_idx) == 0) return(NULL)
    gene_idx <- gene_idx[order(a$start[gene_idx])]
    ge <- a[gene_idx, ]

    in_call <- ge$start <= cnv$End & ge$end >= cnv$Start
    if (!any(in_call)) return(NULL)
    in_idx <- which(in_call)
    row_in_gene <- seq_along(gene_idx)
    exon_type <- case_when(
      in_call                             ~ "in_call",
      row_in_gene == (min(in_idx) - 1L)   ~ "adjacent",
      row_in_gene == (max(in_idx) + 1L)   ~ "adjacent",
      TRUE                                ~ "other"
    )
    target_norm <- norm[gene_idx, target]
    ref_medians <- apply(norm[gene_idx, refs, drop = FALSE], 1, median, na.rm = TRUE)

    tibble(
      call_id     = cnv$call_id,
      Sample      = cnv$Sample,
      run         = cnv$run,
      caller      = cnv$caller,
      Gene        = cnv$Gene,
      CNV_Type    = cnv$CNV_Type,
      Chromosome  = std_chr(cnv$Chromosome),
      CNV_Start   = cnv$Start,
      CNV_End     = cnv$End,
      exon_start  = ge$start,
      exon_end    = ge$end,
      exon_type   = exon_type,
      read_count  = counts$raw[gene_idx, target],
      read_ratio  = round(target_norm / ref_medians, 2)
    )
  })
  bind_rows(out)
}

# RR criteria (corrected version: RR2 checks exons OUTSIDE the call, same gene)
summarise_rr <- function(all_rr_exons, rr_tol = 0.15, rr_del_max = 1.35,
                         rr_dup_min = 0.65, rr_max_spread = 0.5) {
  if (nrow(all_rr_exons) == 0) {
    return(tibble(call_id = character(), n_in_call_exons = integer(),
                  n_out_call_exons = integer(), rr1_pass = logical(),
                  rr2_pass = logical(), rr_spread = numeric(),
                  rr3_pass = logical(), rr_pass_all = logical()))
  }
  safe_spread <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_real_ else max(x) - min(x)
  }
  all_rr_exons %>%
    mutate(
      in_call     = exon_type == "in_call",
      expected_rr = if_else(CNV_Type == "deletion", 0.5, 1.5)
    ) %>%
    group_by(call_id) %>%
    summarise(
      n_in_call_exons  = sum(in_call),
      n_out_call_exons = sum(!in_call),

      # RR1: at least one in-call exon close to the expected ratio
      rr1_pass = any(in_call & abs(read_ratio - expected_rr) <= rr_tol, na.rm = TRUE),

      # RR2: no exon OUTSIDE the call (same gene) with a ratio suggesting the opposite dosage
      rr2_pass = case_when(
        CNV_Type[1L] == "deletion"    ~ !any(!in_call & read_ratio >= rr_del_max, na.rm = TRUE),
        CNV_Type[1L] == "duplication" ~ !any(!in_call & read_ratio <= rr_dup_min, na.rm = TRUE),
        TRUE ~ NA
      ),

      # RR3: spread across in-call exons only
      rr_spread = safe_spread(read_ratio[in_call]),
      rr3_pass  = if_else(n_in_call_exons > 1L, replace_na(rr_spread <= rr_max_spread, FALSE), TRUE),

      rr_pass_all = rr1_pass & rr2_pass & rr3_pass,
      .groups = "drop"
    )
}

# ============================================================
# VAF: WILSON SCORE 95% CI
# ============================================================

wilson_ci <- function(k, n, conf = 0.95) {
  z      <- qnorm(1 - (1 - conf) / 2)
  p_hat  <- k / n
  denom  <- 1 + z^2 / n
  centre <- (p_hat + z^2 / (2 * n)) / denom
  margin <- (z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))) / denom
  list(lower = centre - margin, upper = centre + margin)
}

# Vectorised: does the CI of each VAF include any of `targets`?
ci_includes <- function(vaf, dp, targets) {
  ci <- wilson_ci(round(vaf * dp), dp)
  Reduce(`|`, lapply(targets, function(p) p >= ci$lower & p <= ci$upper))
}

# Duplication:
#   neutral    = CI includes 0.5 AND (1/3 or 2/3)
#   concordant = CI includes 1/3 or 2/3 (not 0.5)
#   discordant = CI excludes 1/3, 2/3
annotate_vaf_status_dup <- function(df) {
  if (nrow(df) == 0) return(df %>% mutate(vaf_status = character()))
  df %>%
    mutate(
      .het = ci_includes(vaf, dp, 0.5),
      .dup = ci_includes(vaf, dp, c(1/3, 2/3)),
      vaf_status = case_when(
        .het & .dup ~ "neutral",
        .dup        ~ "concordant",
        TRUE        ~ "discordant"
      )
    ) %>%
    dplyr::select(-.het, -.dup)
}

# Deletion (flowchart version):
#   discordant = CI overlaps 0.5 (retained heterozygosity)
#   neutral    = otherwise
annotate_vaf_status_del <- function(df) {
  if (nrow(df) == 0) return(df %>% mutate(vaf_status = character()))
  df %>%
    mutate(vaf_status = if_else(ci_includes(vaf, dp, 0.5), "discordant", "neutral"))
}

annotate_vaf_status <- function(df, cnv_type) {
  if (tolower(cnv_type) == "deletion") annotate_vaf_status_del(df) else annotate_vaf_status_dup(df)
}

add_wilson_bounds <- function(df) {
  if (nrow(df) == 0) return(df %>% mutate(ci_lower = numeric(), ci_upper = numeric()))
  ci <- wilson_ci(round(df$vaf * df$dp), df$dp)
  df %>% mutate(ci_lower = ci$lower, ci_upper = ci$upper)
}

# ============================================================
# VCF LOOKUP / SNV EXTRACTION (Rsamtools::scanTabix; no VariantAnnotation)
# ============================================================

build_vcf_lookup <- function(vcf_dir, id_regex = "DNS\\d+", pattern = "\\.sorted\\.vcf\\.gz$") {
  files <- list.files(vcf_dir, pattern = pattern, full.names = TRUE)
  files <- files[file.exists(paste0(files, ".tbi"))]
  tibble(
    vcf_path = files,
    stem     = basename(files) %>% str_remove("\\.vcf\\.gz$") %>% str_remove("\\.sorted$"),
    id       = str_extract(basename(files), id_regex)
  )
}

find_vcf <- function(sample_std, lookup, id_regex = "DNS\\d+") {
  hit <- lookup$vcf_path[lookup$stem == sample_std]
  if (length(hit) == 1) return(hit)
  id  <- str_extract(sample_std, id_regex)
  hit <- lookup$vcf_path[!is.na(lookup$id) & lookup$id == id]
  if (length(hit) == 1) return(hit)
  NA_character_
}

info_field <- function(info, key) {
  m <- str_match(info, paste0("(?:^|;)", key, "=([^;]*)"))[, 2]
  m
}

# SNVs in region with DP >= dp_min. VAF from INFO/AF (falls back to FORMAT/AD).
# DP from FORMAT/DP (falls back to INFO/DP).
extract_snvs_in_region <- function(vcf_path, chrom, start, end, dp_min = 30) {
  empty <- tibble(pos = integer(), ref = character(), alt = character(),
                  vaf = numeric(), dp = numeric())
  if (is.na(vcf_path) || !file.exists(vcf_path)) return(empty)
  tryCatch({
    tbx  <- TabixFile(vcf_path)
    seqs <- seqnamesTabix(tbx)
    chr  <- std_chr(chrom)
    seqn <- if (paste0("chr", chr) %in% seqs) paste0("chr", chr) else if (chr %in% seqs) chr else NA
    if (is.na(seqn)) return(empty)
    lines <- scanTabix(tbx, param = GRanges(seqn, IRanges(start, end)))[[1]]
    if (length(lines) == 0) return(empty)

    f      <- str_split_fixed(lines, "\t", 10)
    ref    <- f[, 4]
    alt    <- f[, 5]
    info   <- f[, 8]
    fmt    <- str_split(f[, 9], ":")
    smp    <- str_split(f[, 10], ":")
    fmt_get <- function(key) mapply(function(k, v) {
      j <- match(key, k); if (is.na(j) || j > length(v)) NA_character_ else v[j]
    }, fmt, smp, USE.NAMES = FALSE)

    type   <- info_field(info, "TYPE")
    is_snv <- if (all(is.na(type))) {
      nchar(ref) == 1 & grepl("^[ACGTN](,[ACGTN])*$", alt)
    } else {
      !is.na(type) & type == "SNV"
    }

    af <- suppressWarnings(as.numeric(sub(",.*", "", info_field(info, "AF"))))
    if (all(is.na(af))) {
      ad  <- str_split(fmt_get("AD"), ",")
      af  <- vapply(ad, function(x) {
        x <- suppressWarnings(as.numeric(x)); if (length(x) < 2 || sum(x) == 0) NA_real_ else sum(x[-1]) / sum(x)
      }, numeric(1))
    }
    dp <- suppressWarnings(as.numeric(fmt_get("DP")))
    if (all(is.na(dp))) dp <- suppressWarnings(as.numeric(info_field(info, "DP")))

    tibble(pos = as.integer(f[, 2]), ref = ref, alt = alt, vaf = af, dp = dp) %>%
      filter(is_snv, !is.na(vaf), !is.na(dp), dp >= dp_min)
  }, error = function(e) {
    warning(sprintf("VCF read error %s at %s:%d-%d - %s", vcf_path, chrom, start, end, conditionMessage(e)))
    empty
  })
}

# Informative SNVs for a set of calls, annotated per the call's CNV type.
# Duplications: 0 < VAF < vaf_max ; Deletions: 0 < VAF < vaf_max
extract_vafs_for_calls <- function(calls_df, vcf_lookup, vaf_max = 0.95, dp_min = 30,
                                   id_regex = "DNS\\d+") {
  out <- lapply(seq_len(nrow(calls_df)), function(i) {
    cnv <- calls_df[i, ]
    vcf <- find_vcf(cnv$sample_std, vcf_lookup, id_regex)
    snv <- extract_snvs_in_region(vcf, cnv$Chromosome, cnv$Start, cnv$End, dp_min)
    if (nrow(snv) == 0) return(NULL)
    snv %>%
      filter(vaf > 0, vaf < vaf_max) %>%
      annotate_vaf_status(cnv$CNV_Type) %>%
      mutate(call_id = cnv$call_id, Sample = cnv$Sample, run = cnv$run,
             caller = cnv$caller, Gene = cnv$Gene, CNV_Type = cnv$CNV_Type,
             Chromosome = std_chr(cnv$Chromosome), CNV_Start = cnv$Start,
             CNV_End = cnv$End, vcf_path = vcf, .before = 1)
  })
  res <- bind_rows(out)
  if (ncol(res) == 0) {
    res <- tibble(call_id = character(), Sample = character(), run = character(),
                  caller = character(), Gene = character(), CNV_Type = character(),
                  Chromosome = character(), CNV_Start = integer(), CNV_End = integer(),
                  vcf_path = character(), pos = integer(), ref = character(),
                  alt = character(), vaf = numeric(), dp = numeric(),
                  vaf_status = character())
  }
  res
}
