#!/usr/bin/env Rscript
# ============================================================
# STAGE 2: PRELIMINARY FILTERING
#   DECoN:    BF >= bf_min (all calls)
#             AND mean read count >= min_rc (single-exon calls only)
#   clearCNV: mean read count >= min_rc (single-exon calls only)
# Single-exon = the call overlaps exactly one framework-BED exon of its gene.
# Mean read count = mean raw count over all BED exons overlapping the call.
# Output keeps every raw call, with pass/fail columns.
# ============================================================

suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--run"),
  make_option("--decon",    help = "DECoN calls_all.txt"),
  make_option("--clearcnv", help = "clearCNV cnv_calls.tsv"),
  make_option("--counts",   help = "counts RData"),
  make_option("--outliers", help = "<run>_outliers.tsv"),
  make_option("--bed",      help = "framework BED"),
  make_option("--bf_min",   type = "double", default = 10),
  make_option("--min_rc",   type = "double", default = 100),
  make_option("--utils")
)))
source(opt$utils)
run <- opt$run

# ── Standardise ──────────────────────────────────────────────────────────────

decon_raw <- read_tsv_safe(opt$decon)
decon <- if (nrow(decon_raw) == 0) tibble() else decon_raw %>%
  transmute(
    caller     = "DECoN",
    Sample     = std_sample(Sample),
    Chromosome = std_chr(Chromosome),
    Start      = as.integer(Start),
    End        = as.integer(End),
    Gene       = as.character(Gene),
    CNV_Type   = as.character(CNV.type),
    Quality    = as.numeric(BF),
    caller_n_exons = as.integer(N.exons),
    caller_ratio   = as.numeric(Reads.ratio)
  )

clear_raw <- read_tsv_safe(opt$clearcnv)
clear <- if (nrow(clear_raw) == 0) tibble() else clear_raw %>%
  transmute(
    caller     = "clearCNV",
    Sample     = std_sample(sample),
    Chromosome = std_chr(chr),
    Start      = as.integer(start),
    End        = as.integer(end),
    Gene       = as.character(gene),
    CNV_Type   = recode(as.character(aberration), "DEL" = "deletion", "DUP" = "duplication"),
    Quality    = as.numeric(score),
    caller_n_exons = as.integer(size),
    caller_ratio   = as.numeric(ratio)
  )

template <- tibble(caller = character(), Sample = character(), Chromosome = character(),
                   Start = integer(), End = integer(), Gene = character(),
                   CNV_Type = character(), Quality = numeric(),
                   caller_n_exons = integer(), caller_ratio = numeric())

calls <- bind_rows(template, decon, clear) %>%
  mutate(
    run               = run,
    sample_std        = Sample,
    sample_id_trimmed = str_extract(Sample, "DNS\\d+"),
    call_id           = make_call_id(run, caller, Sample, Chromosome, Start, End, CNV_Type),
    .before = 1
  )

# ── Features ────────────────────────────────────────────────────────────────

bed      <- read_framework_bed(opt$bed)
counts   <- load_counts(opt$counts)
outliers <- read_tsv_safe(opt$outliers)
outlier_samples <- outliers$sample[outliers$is_outlier %in% TRUE]

calls <- count_exons_in_call(calls, bed)
calls$mean_read_count <- call_mean_read_count(calls, counts)

# ── Filters ─────────────────────────────────────────────────────────────────

calls <- calls %>%
  mutate(
    is_single_exon = n_exons_in_call == 1L,
    bf_pass = if_else(caller == "DECoN", replace_na(Quality >= opt$bf_min, FALSE), TRUE),
    rc_pass = if_else(is_single_exon, replace_na(mean_read_count >= opt$min_rc, FALSE), TRUE),
    from_outlier = sample_std %in% outlier_samples,
    prelim_pass = bf_pass & rc_pass & !from_outlier,
    prelim_fail_reason = case_when(
      from_outlier         ~ "outlier_sample",
      !bf_pass & !rc_pass  ~ "low_BF;low_read_count_single_exon",
      !bf_pass             ~ "low_BF",
      !rc_pass             ~ "low_read_count_single_exon",
      TRUE                 ~ NA_character_
    )
  )

write_tsv_safe(calls, sprintf("%s_stage2_calls.tsv", run))
message(sprintf("%s: %d raw calls (DECoN %d, clearCNV %d); %d pass preliminary filter",
                run, nrow(calls), sum(calls$caller == "DECoN"),
                sum(calls$caller == "clearCNV"), sum(calls$prelim_pass)))
