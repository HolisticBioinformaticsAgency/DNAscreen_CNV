#!/usr/bin/env Rscript
# ============================================================
# STAGE 4: ASSESSMENT OF READ RATIO AND VARIANT ALLELE FRACTION
# Decision applied to single-caller calls only; RR/VAF features are
# also computed for dual-caller calls (for tracking and plots).
#
# RR (reference = top-N retained samples by average percent difference):
#   RR1  >= 1 in-call exon within +/- rr_tol of expected (0.5 del / 1.5 dup)
#   RR2  no exon OUTSIDE the call (same gene) with RR >= rr_del_max (del)
#        or RR <= rr_dup_min (dup)
#   RR3  max - min RR across in-call exons <= rr_max_spread (multi-exon)
#
# VAF (SNVs, 0 < VAF < vaf_max, DP >= dp_min, Wilson 95% CI):
#   Duplication  no SNV:  single-exon -> pass ; multi-exon -> fail
#                SNVs:    pass if prop concordant (CI overlaps 1/3 or 2/3,
#                         not 0.5) > dup_min_prop AND n discordant <= dup_max_disc
#   Deletion     multi-exon -> fail (VAFs cannot support it)
#                single-exon: pass if no SNV CI overlaps 0.5 (or no SNVs)
# ============================================================

suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--run"),
  make_option("--stage3"),
  make_option("--counts"),
  make_option("--outliers"),
  make_option("--vcf_dir"),
  make_option("--vcf_pattern",   default = "\\.sorted\\.vcf\\.gz$"),
  make_option("--id_regex",      default = "DNS\\d+"),
  make_option("--n_ref",         type = "integer", default = 20L),
  make_option("--rr_tol",        type = "double",  default = 0.15),
  make_option("--rr_del_max",    type = "double",  default = 1.35),
  make_option("--rr_dup_min",    type = "double",  default = 0.65),
  make_option("--rr_max_spread", type = "double",  default = 0.5),
  make_option("--vaf_max",       type = "double",  default = 0.95),
  make_option("--dp_min",        type = "double",  default = 30),
  make_option("--dup_min_prop",  type = "double",  default = 0.5),
  make_option("--dup_max_disc",  type = "integer", default = 2L),
  make_option("--utils")
)))
source(opt$utils)
run <- opt$run

calls    <- read_tsv_safe(opt$stage3)
outliers <- read_tsv_safe(opt$outliers)
counts   <- load_counts(opt$counts)

retained <- outliers$sample[!(outliers$is_outlier %in% TRUE)]
retained <- intersect(retained, colnames(counts$raw))
norm     <- cpt_normalise(counts$raw[, retained, drop = FALSE])
apd      <- compute_apd(norm)

if (nrow(calls) == 0) {
  write_tsv_safe(tibble(), sprintf("%s_stage4_calls.tsv", run))
  write_tsv_safe(tibble(), sprintf("%s_rr_per_exon.tsv", run))
  write_tsv_safe(tibble(), sprintf("%s_vaf_per_snv.tsv", run))
  quit(save = "no")
}

# ── Read ratios ─────────────────────────────────────────────────────────────
rr_exons <- exon_read_ratios(calls, counts, norm, apd, n_ref = opt$n_ref)
rr_summary <- summarise_rr(rr_exons, opt$rr_tol, opt$rr_del_max,
                           opt$rr_dup_min, opt$rr_max_spread)

# ── VAFs ────────────────────────────────────────────────────────────────────
vcf_lookup <- build_vcf_lookup(opt$vcf_dir, opt$id_regex, opt$vcf_pattern)
vafs <- extract_vafs_for_calls(calls, vcf_lookup, opt$vaf_max, opt$dp_min, opt$id_regex)
vafs <- add_wilson_bounds(vafs)

vaf_summary <- vafs %>%
  group_by(call_id) %>%
  summarise(
    n_informative_vafs = n(),
    n_concordant       = sum(vaf_status == "concordant", na.rm = TRUE),
    n_neutral          = sum(vaf_status == "neutral",    na.rm = TRUE),
    n_discordant       = sum(vaf_status == "discordant", na.rm = TRUE),
    .groups = "drop"
  )

calls$vcf_found <- vapply(calls$sample_std, function(s)
  !is.na(find_vcf(s, vcf_lookup, opt$id_regex)), logical(1))

# ── Decision rules ──────────────────────────────────────────────────────────
stage4 <- calls %>%
  left_join(rr_summary,  by = "call_id") %>%
  left_join(vaf_summary, by = "call_id") %>%
  mutate(
    n_informative_vafs = replace_na(n_informative_vafs, 0L),
    n_concordant       = replace_na(n_concordant, 0L),
    n_neutral          = replace_na(n_neutral, 0L),
    n_discordant       = replace_na(n_discordant, 0L),
    prop_concordant    = if_else(n_informative_vafs > 0,
                                 n_concordant / n_informative_vafs, NA_real_),
    is_single_exon     = replace_na(n_exons_in_call == 1L, FALSE),
    rr_pass_all        = replace_na(rr_pass_all, FALSE),

    vaf_reason = case_when(
      # ── DELETION ──
      CNV_Type == "deletion" & !is_single_exon            ~ "fail:multi_exon_deletion",
      CNV_Type == "deletion" & n_informative_vafs == 0L   ~ "pass:no_informative_snv",
      CNV_Type == "deletion" & n_discordant > 0L          ~ "fail:snv_ci_overlaps_0.5",
      CNV_Type == "deletion"                              ~ "pass:no_snv_ci_overlaps_0.5",
      # ── DUPLICATION ──
      CNV_Type == "duplication" & n_informative_vafs == 0L & !is_single_exon ~ "fail:no_informative_snv_multi_exon",
      CNV_Type == "duplication" & n_informative_vafs == 0L                   ~ "pass:no_informative_snv",
      CNV_Type == "duplication" & prop_concordant > opt$dup_min_prop &
        n_discordant <= opt$dup_max_disc                                      ~ "pass:vafs_support_duplication",
      CNV_Type == "duplication"                                               ~ "fail:vafs_do_not_support_duplication",
      TRUE ~ "fail:unknown_cnv_type"
    ),
    vaf_pass = startsWith(vaf_reason, "pass"),

    rr_reason = case_when(
      is.na(n_in_call_exons) ~ "fail:no_gene_exons_in_counts",
      rr_pass_all            ~ "pass",
      TRUE ~ paste0("fail:", paste(
        if_else(rr1_pass %in% TRUE, "", "RR1"),
        if_else(rr2_pass %in% TRUE, "", "RR2"),
        if_else(rr3_pass %in% TRUE, "", "RR3"), sep = ","
      ) %>% str_replace_all("^,+|,+$", "") %>% str_replace_all(",+", ","))
    ),

    stage4_status = case_when(
      role == "anchor"  ~ "not_assessed:dual_caller",
      role == "merged"  ~ "not_assessed:dual_caller",
      !rr_pass_all      ~ "rejected_rr",
      !vaf_pass         ~ "rejected_vaf",
      TRUE              ~ "passed_rr_vaf"
    ),

    final_status = case_when(
      role == "anchor"                  ~ "high_confidence_dual_caller",
      role == "merged"                  ~ "merged_into_dual_caller",
      stage4_status == "passed_rr_vaf"  ~ "high_confidence_rr_vaf",
      TRUE                              ~ stage4_status
    )
  )

write_tsv_safe(stage4,   sprintf("%s_stage4_calls.tsv", run))
write_tsv_safe(rr_exons, sprintf("%s_rr_per_exon.tsv", run))
write_tsv_safe(vafs,     sprintf("%s_vaf_per_snv.tsv", run))

single <- stage4 %>% filter(role == "single_caller")
message(sprintf("%s: single-caller calls %d -> passed %d, rejected RR %d, rejected VAF %d",
                run, nrow(single), sum(single$stage4_status == "passed_rr_vaf"),
                sum(single$stage4_status == "rejected_rr"),
                sum(single$stage4_status == "rejected_vaf")))
