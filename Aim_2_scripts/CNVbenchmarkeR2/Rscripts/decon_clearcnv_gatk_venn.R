library(dplyr)
library(GenomicRanges)
library(ggVennDiagram)

decon_calls <- harmonised_list_with_acgh$decon$data %>% filter(acgh_status != "tbc") 
clearcnv_calls <- harmonised_list_with_acgh$clearCNV$data %>% filter(acgh_status != "tbc") 
gatk_calls <- harmonised_list_with_acgh$gatk$data %>% filter(acgh_status != "tbc")

decon_calls <- decon_calls %>%
  filter(BF >= 10)

gatk_calls <- gatk_calls %>%
  filter(gt_QS >= 160)

# ---- Convert calls dataframe to GRanges with Sample info ----
df_to_gr <- function(df) {
  GRanges(
    seqnames = df$Chromosome,
    ranges = IRanges(start = df$Start, end = df$End),
    Sample = df$Sample,
    acgh_status = df$acgh_status
  )
}

decon_gr <- df_to_gr(decon_calls)
clearcnv_gr <- df_to_gr(clearcnv_calls)
gatk_gr <- df_to_gr(gatk_calls)

# ---- Function to find overlapping calls between two GRanges with Sample match ----
get_overlaps <- function(gr1, gr2) {
  hits <- findOverlaps(gr1, gr2)
  hits <- hits[gr1$Sample[queryHits(hits)] == gr2$Sample[subjectHits(hits)]]
  df1_hits <- gr1[queryHits(hits)]
  df2_hits <- gr2[subjectHits(hits)]
  data.frame(
    Sample = df1_hits$Sample,
    acgh_status1 = df1_hits$acgh_status,
    acgh_status2 = df2_hits$acgh_status
  )
}

# ---- Pairwise overlaps ----
decon_clearcnv <- get_overlaps(decon_gr, clearcnv_gr)
decon_gatk <- get_overlaps(decon_gr, gatk_gr)
clearcnv_gatk <- get_overlaps(clearcnv_gr, gatk_gr)

# ---- Find trio overlaps ----
dc_hits <- findOverlaps(decon_gr, clearcnv_gr)
dc_hits <- dc_hits[decon_gr$Sample[queryHits(dc_hits)] == clearcnv_gr$Sample[subjectHits(dc_hits)]]

trio_list <- lapply(seq_along(queryHits(dc_hits)), function(i) {
  decon_idx <- queryHits(dc_hits)[i]
  clearcnv_idx <- subjectHits(dc_hits)[i]
  
  g_hits <- findOverlaps(decon_gr[decon_idx], gatk_gr)
  g_hits <- g_hits[decon_gr$Sample[decon_idx] == gatk_gr$Sample[subjectHits(g_hits)]]
  
  if (length(g_hits) > 0) {
    data.frame(
      Sample = decon_gr$Sample[decon_idx],
      acgh_status_decon = decon_gr$acgh_status[decon_idx],
      acgh_status_clearcnv = clearcnv_gr$acgh_status[clearcnv_idx],
      acgh_status_gatk = gatk_gr$acgh_status[subjectHits(g_hits)]
    )
  } else {
    NULL
  }
})

trio_hits <- do.call(rbind, trio_list)

# ---- Find tool-specific calls (not shared with others) ----
get_unique_calls <- function(gr1, gr2, gr3) {
  hits2 <- findOverlaps(gr1, gr2)
  hits2 <- hits2[gr1$Sample[queryHits(hits2)] == gr2$Sample[subjectHits(hits2)]]
  
  hits3 <- findOverlaps(gr1, gr3)
  hits3 <- hits3[gr1$Sample[queryHits(hits3)] == gr3$Sample[subjectHits(hits3)]]
  
  unique_idx <- setdiff(seq_along(gr1), unique(c(queryHits(hits2), queryHits(hits3))))
  gr1[unique_idx]
}

decon_only <- get_unique_calls(decon_gr, clearcnv_gr, gatk_gr)
clearcnv_only <- get_unique_calls(clearcnv_gr, decon_gr, gatk_gr)
gatk_only <- get_unique_calls(gatk_gr, decon_gr, clearcnv_gr)

# ---- Compute pairwise overlaps excluding trio (Decon & GATK only, etc.) ----
decon_gatk_only <- decon_gatk %>%
  anti_join(trio_hits %>% select(Sample, acgh_status_decon, acgh_status_gatk),
            by = c("Sample" = "Sample",
                   "acgh_status1" = "acgh_status_decon",
                   "acgh_status2" = "acgh_status_gatk"))

decon_clearcnv_only <- decon_clearcnv %>%
  anti_join(trio_hits %>% select(Sample, acgh_status_decon, acgh_status_clearcnv),
            by = c("Sample" = "Sample",
                   "acgh_status1" = "acgh_status_decon",
                   "acgh_status2" = "acgh_status_clearcnv"))

clearcnv_gatk_only <- clearcnv_gatk %>%
  anti_join(trio_hits %>% select(Sample, acgh_status_clearcnv, acgh_status_gatk),
            by = c("Sample" = "Sample",
                   "acgh_status1" = "acgh_status_clearcnv",
                   "acgh_status2" = "acgh_status_gatk"))

# ---- Inspect results ----
decon_gatk_only
decon_clearcnv_only
clearcnv_gatk_only
trio_hits
decon_only
clearcnv_only
gatk_only
