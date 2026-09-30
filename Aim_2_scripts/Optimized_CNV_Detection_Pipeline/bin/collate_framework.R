#!/usr/bin/env Rscript
# ============================================================
# COLLATE: track BAMs and CNV calls through the framework, across runs
# Reads every per-run table staged in the working directory and writes:
#   01_outliers/ 02_raw_calls/ 03_prelim_filter/ 04_prioritisation/
#   05_rr_vaf/ 06_final/ call_ledger.tsv sample_ledger.tsv
#   framework_summary.tsv framework_summary_by_run.tsv run_status.tsv
# ============================================================

suppressPackageStartupMessages(library(optparse))
opt <- parse_args(OptionParser(option_list = list(make_option("--utils"))))
source(opt$utils)

read_all <- function(pattern) {
  files <- list.files(".", pattern = pattern, full.names = TRUE)
  dfs <- lapply(files, function(f) {
    if (file.size(f) == 0) return(NULL)
    suppressMessages(read_tsv(f, col_types = cols(.default = col_character()), na = c("", "NA")))
  })
  out <- bind_rows(dfs)
  if (nrow(out) == 0) return(out)
  suppressMessages(type_convert(out, guess_integer = TRUE, na = c("", "NA")))
}

run_order <- function(x) order(suppressWarnings(as.numeric(str_extract(x, "\\d+"))), x)
sort_runs <- function(df) if (nrow(df) == 0) df else df[run_order(df$run), ]

outliers <- read_all("_outliers\\.tsv$")   %>% sort_runs()
stage2   <- read_all("_stage2_calls\\.tsv$") %>% sort_runs()
stage4   <- read_all("_stage4_calls\\.tsv$") %>% sort_runs()
rr_exons <- read_all("_rr_per_exon\\.tsv$")  %>% sort_runs()
vafs     <- read_all("_vaf_per_snv\\.tsv$")  %>% sort_runs()
run_stat <- read_all("_caller_status\\.tsv$") %>% sort_runs()

if (nrow(stage2) == 0) {
  stage2 <- tibble(call_id = character(), run = character(), caller = character(),
                   sample_std = character(), prelim_pass = logical(),
                   prelim_fail_reason = character())
}

core_cols <- c("call_id", "run", "caller", "Sample", "sample_id_trimmed", "Chromosome",
               "Start", "End", "Gene", "CNV_Type", "Quality", "n_exons_in_call",
               "mean_read_count")

# ── 01 outliers ─────────────────────────────────────────────────────────────
write_tsv_safe(outliers %>% filter(is_outlier %in% TRUE), "01_outliers/outlier_bams.tsv")
write_tsv_safe(outliers %>% filter(!(is_outlier %in% TRUE)), "01_outliers/retained_bams.tsv")

# ── 02 raw calls ────────────────────────────────────────────────────────────
write_tsv_safe(stage2 %>% filter(caller == "DECoN"),    "02_raw_calls/decon_calls.tsv")
write_tsv_safe(stage2 %>% filter(caller == "clearCNV"), "02_raw_calls/clearcnv_calls.tsv")

# ── 03 preliminary filter ───────────────────────────────────────────────────
write_tsv_safe(stage2 %>% filter(prelim_pass %in% TRUE),    "03_prelim_filter/passed.tsv")
write_tsv_safe(stage2 %>% filter(!(prelim_pass %in% TRUE)), "03_prelim_filter/failed.tsv")

# ── 04 prioritisation ───────────────────────────────────────────────────────
if (nrow(stage4) == 0) {
  stage4 <- tibble(call_id = character(), role = character(), priority = character(),
                   stage4_status = character(), final_status = character(),
                   partner_call_id = character(), reciprocal_overlap = numeric(),
                   rr_reason = character(), vaf_reason = character())
}
anchors <- stage4 %>% filter(role == "anchor")
dual <- anchors %>%
  left_join(stage4 %>% filter(role == "merged") %>%
              dplyr::select(partner_call_id = call_id, clearcnv_Start = Start,
                            clearcnv_End = End, clearcnv_Quality = Quality),
            by = "partner_call_id")
write_tsv_safe(dual, "04_prioritisation/dual_caller_high_confidence.tsv")
write_tsv_safe(stage4 %>% filter(role == "merged"), "04_prioritisation/dual_caller_clearcnv_partners.tsv")
single <- stage4 %>% filter(role == "single_caller")
write_tsv_safe(single, "04_prioritisation/single_caller_manual_inspection.tsv")

# ── 05 RR + VAF ─────────────────────────────────────────────────────────────
write_tsv_safe(single %>% filter(stage4_status == "passed_rr_vaf"), "05_rr_vaf/passed.tsv")
write_tsv_safe(single %>% filter(stage4_status == "rejected_rr"),   "05_rr_vaf/failed_rr.tsv")
write_tsv_safe(single %>% filter(stage4_status == "rejected_vaf"),  "05_rr_vaf/failed_vaf.tsv")
write_tsv_safe(rr_exons, "05_rr_vaf/rr_per_exon.tsv")
write_tsv_safe(vafs,     "05_rr_vaf/vaf_per_snv.tsv")

# ── Call ledger ─────────────────────────────────────────────────────────────
s4_cols <- intersect(c("call_id", "priority", "role", "partner_call_id", "reciprocal_overlap",
                       "n_in_call_exons", "n_out_call_exons", "rr1_pass", "rr2_pass",
                       "rr_spread", "rr3_pass", "rr_pass_all", "rr_reason",
                       "vcf_found", "n_informative_vafs", "n_concordant", "n_neutral",
                       "n_discordant", "prop_concordant", "vaf_pass", "vaf_reason",
                       "stage4_status", "final_status"), names(stage4))

ledger <- stage2 %>%
  left_join(stage4 %>% dplyr::select(all_of(s4_cols)), by = "call_id") %>%
  mutate(
    final_status = if_else(prelim_pass %in% TRUE, final_status, "rejected_prelim_filter"),
    stage_reached = case_when(
      !(prelim_pass %in% TRUE)                 ~ "2_prelim_filter",
      role %in% c("anchor", "merged")          ~ "3_prioritisation",
      TRUE                                     ~ "4_rr_vaf"
    ),
    fail_reason = case_when(
      final_status == "rejected_prelim_filter" ~ prelim_fail_reason,
      final_status == "rejected_rr"            ~ rr_reason,
      final_status == "rejected_vaf"           ~ vaf_reason,
      TRUE                                     ~ NA_character_
    ),
    is_high_confidence = final_status %in% c("high_confidence_dual_caller",
                                             "high_confidence_rr_vaf",
                                             "merged_into_dual_caller")
  )
write_tsv_safe(ledger, "call_ledger.tsv")

# ── 06 final ────────────────────────────────────────────────────────────────
final_hc <- ledger %>%
  filter(final_status %in% c("high_confidence_dual_caller", "high_confidence_rr_vaf")) %>%
  mutate(called_by = if_else(final_status == "high_confidence_dual_caller",
                             "DECoN+clearCNV", caller))
write_tsv_safe(final_hc, "06_final/high_confidence_calls.tsv")
write_tsv_safe(ledger %>% filter(str_starts(final_status, "rejected")), "06_final/rejected_calls.tsv")

# ── Sample ledger ───────────────────────────────────────────────────────────
per_sample <- ledger %>%
  group_by(run, sample = sample_std) %>%
  summarise(
    n_raw_decon      = sum(caller == "DECoN"),
    n_raw_clearcnv   = sum(caller == "clearCNV"),
    n_prelim_pass    = sum(prelim_pass %in% TRUE),
    n_dual_caller    = sum(role %in% "anchor"),
    n_single_caller  = sum(role %in% "single_caller"),
    n_rejected_rr    = sum(final_status == "rejected_rr"),
    n_rejected_vaf   = sum(final_status == "rejected_vaf"),
    n_high_confidence = sum(final_status %in% c("high_confidence_dual_caller", "high_confidence_rr_vaf")),
    high_confidence_calls = paste(call_id[final_status %in% c("high_confidence_dual_caller",
                                                               "high_confidence_rr_vaf")], collapse = ";"),
    .groups = "drop"
  )
sample_ledger <- outliers %>%
  dplyr::select(run, sample, bam, total_reads, max_correlation, cluster, is_outlier, outlier_reason) %>%
  left_join(per_sample, by = c("run", "sample")) %>%
  mutate(across(starts_with("n_"), ~ replace_na(.x, 0L)),
         high_confidence_calls = na_if(high_confidence_calls, ""))
write_tsv_safe(sample_ledger, "sample_ledger.tsv")

# ── Summary counts (the numbers on the framework figure) ────────────────────
summarise_counts <- function(led, outl) {
  bind_rows(
    tibble(stage = "1_outliers", caller = "all", category = "bams_total",    n = nrow(outl)),
    tibble(stage = "1_outliers", caller = "all", category = "bams_outlier",  n = sum(outl$is_outlier %in% TRUE)),
    tibble(stage = "1_outliers", caller = "all", category = "bams_retained", n = sum(!(outl$is_outlier %in% TRUE))),
    led %>% dplyr::count(caller) %>% mutate(stage = "2_identify_cnvs", category = "raw_calls"),
    led %>% filter(prelim_pass %in% TRUE) %>% dplyr::count(caller) %>%
      mutate(stage = "3_prelim_filter", category = "passed"),
    led %>% filter(!(prelim_pass %in% TRUE)) %>% dplyr::count(caller, category = prelim_fail_reason) %>%
      mutate(stage = "3_prelim_filter", category = paste0("failed:", category)),
    tibble(stage = "4_prioritisation", caller = "DECoN+clearCNV", category = "dual_caller",
           n = sum(led$role %in% "anchor")),
    led %>% filter(role %in% "single_caller") %>% dplyr::count(caller) %>%
      mutate(stage = "4_prioritisation", category = "single_caller"),
    led %>% filter(role %in% "single_caller") %>% dplyr::count(caller, category = final_status) %>%
      mutate(stage = "5_rr_vaf"),
    tibble(stage = "6_final", caller = "all", category = "high_confidence",
           n = sum(led$final_status %in% c("high_confidence_dual_caller", "high_confidence_rr_vaf")))
  ) %>% dplyr::select(stage, caller, category, n)
}

write_tsv_safe(summarise_counts(ledger, outliers), "framework_summary.tsv")

by_run <- bind_rows(lapply(unique(c(outliers$run, ledger$run)), function(r) {
  summarise_counts(ledger %>% filter(run == r), outliers %>% filter(run == r)) %>%
    mutate(run = r, .before = 1)
}))
write_tsv_safe(sort_runs(by_run), "framework_summary_by_run.tsv")
write_tsv_safe(run_stat, "run_status.tsv")

print(as.data.frame(summarise_counts(ledger, outliers)))
