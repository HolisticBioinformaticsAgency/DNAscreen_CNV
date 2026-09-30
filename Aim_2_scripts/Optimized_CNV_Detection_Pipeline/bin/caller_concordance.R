#!/usr/bin/env Rscript
# ============================================================
# STAGE 3: PRIORITISATION MATRIX (caller concordance)
#   Same sample + same CNV type + reciprocal overlap >= threshold
#   -> dual-caller (high confidence, anchored on the DECoN call;
#      the matched clearCNV call is recorded as merged into it)
#   Otherwise -> single-caller (needs RR + VAF assessment)
# Input: stage-2 table (only prelim_pass calls are considered).
# ============================================================

suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--run"),
  make_option("--stage2"),
  make_option("--reciprocal", type = "double", default = 0.5),
  make_option("--utils")
)))
source(opt$utils)
run <- opt$run

calls <- read_tsv_safe(opt$stage2)
if (nrow(calls) > 0) calls <- calls %>% filter(prelim_pass %in% TRUE)

calls <- calls %>%
  mutate(priority = NA_character_, role = NA_character_,
         partner_call_id = NA_character_, partner_quality = NA_real_,
         reciprocal_overlap = NA_real_)

decon <- calls %>% filter(caller == "DECoN")
clear <- calls %>% filter(caller == "clearCNV")

if (nrow(decon) > 0 && nrow(clear) > 0) {
  gr_d <- GRanges(decon$Chromosome, IRanges(decon$Start, decon$End))
  gr_c <- GRanges(clear$Chromosome, IRanges(clear$Start, clear$End))
  hits <- findOverlaps(gr_d, gr_c, ignore.strand = TRUE)
  q <- queryHits(hits); s <- subjectHits(hits)

  keep <- decon$sample_std[q] == clear$sample_std[s] &
    decon$CNV_Type[q] == clear$CNV_Type[s]
  q <- q[keep]; s <- s[keep]

  if (length(q) > 0) {
    inter <- width(pintersect(ranges(gr_d)[q], ranges(gr_c)[s]))
    recip <- pmin(inter / width(gr_d)[q], inter / width(gr_c)[s])
    pairs <- tibble(q = q, s = s, recip = recip) %>% filter(recip >= opt$reciprocal)

    if (nrow(pairs) > 0) {
      # Best clearCNV match per DECoN call
      best_d <- pairs %>% group_by(q) %>%
        slice_max(recip, n = 1, with_ties = FALSE) %>% ungroup()
      decon$priority[best_d$q]           <- "dual_caller"
      decon$role[best_d$q]               <- "anchor"
      decon$partner_call_id[best_d$q]    <- clear$call_id[best_d$s]
      decon$partner_quality[best_d$q]    <- clear$Quality[best_d$s]
      decon$reciprocal_overlap[best_d$q] <- best_d$recip

      # Every matched clearCNV call is merged into its best DECoN match
      best_c <- pairs %>% group_by(s) %>%
        slice_max(recip, n = 1, with_ties = FALSE) %>% ungroup()
      clear$priority[best_c$s]           <- "dual_caller"
      clear$role[best_c$s]               <- "merged"
      clear$partner_call_id[best_c$s]    <- decon$call_id[best_c$q]
      clear$partner_quality[best_c$s]    <- decon$Quality[best_c$q]
      clear$reciprocal_overlap[best_c$s] <- best_c$recip
    }
  }
}

out <- bind_rows(decon, clear) %>%
  mutate(
    priority = coalesce(priority, if_else(caller == "DECoN", "decon_only", "clearcnv_only")),
    role     = coalesce(role, "single_caller")
  )

write_tsv_safe(out, sprintf("%s_stage3_calls.tsv", run))
message(sprintf("%s: %d dual-caller (DECoN anchors), %d DECoN-only, %d clearCNV-only",
                run, sum(out$role == "anchor"), sum(out$priority == "decon_only"),
                sum(out$priority == "clearcnv_only")))
