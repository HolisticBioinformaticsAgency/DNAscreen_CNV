library(dplyr)
library(GenomicRanges)
library(ggVennDiagram)
library(stringr)

# ---- Standardise DeCON ----
standardise_decon <- function(df) {
  df %>%
    transmute(
      Sample,
      sample_id_trimmed = sub(".*-(DNS\\d+)_S\\d+.*", "\\1", Sample),
      run = as.numeric(gsub("run", "", run)),
      Chromosome = as.character(Chromosome),   # keep as string for consistency
      Start = Start,
      End = End,
      Gene = Gene,
      CNV_Type = CNV.type,   # already "deletion"/"duplication"
      Quality = BF
    )
}

# ---- Standardise ClearCNV ----
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


###### only doing for ten concordant runs for now ####
# concordant_runs <- c(31, 14, 25, 8, 19, 4, 21, 20, 13, 16)
# # Function to filter a tool's data by concordant runs
# filter_concordant_runs <- function(df) {
#   df %>%
#     mutate(run_numeric = as.numeric(str_remove(run, "^run"))) %>%
#     filter(run_numeric %in% concordant_runs)
# }

## note clearcnv was run at 0.005 transition probability instead of 0.001 because
## when running with the covered.bed file, we can pick the LDLR dup in DNSxxxxxx up
## at 0.005. So we changed the trans prob to 0.005. But... when running with the
## exon bed, it can't pick it up again... we didn't want to increase the trans prob
## more than necessary so just left it at 0.005.
## UPDATE: we reduced the trans prob to 0.001 using exon bed as clearCNV is picking the uncertain
## DNSxxxxxx call when it is at 0.005...


decon_file <- read.csv("/path/to/project/Aim_2/Optimized_CNV_Detection_Pipeline/output/combined_cnv_calls_decon_exon_bed.tsv", sep="\t")
clearcnv_file <- read.csv("/path/to/project/Aim_2/Optimized_CNV_Detection_Pipeline/output/combined_cnv_calls_clearcnv_exon_bed_001.tsv", sep="\t")


# ---- Standardise and filter ----
perrun_decon_calls <- standardise_decon(decon_file)
perrun_decon_calls$Sample <- sub("\\.hq.sorted.marked*", "", perrun_decon_calls$Sample)

perrun_clearcnv_calls <- standardise_clearcnv(clearcnv_file)
perrun_clearcnv_calls$Sample <- sub("\\.hq.sorted.marked*", "", perrun_clearcnv_calls$Sample)


# Assign aCGH concordance status and non-vscnv status
# Note that for decon and clearcnv, we checked that there are NO cnv calls
# (from the 77 array samples) that are not called by vs-cnv that overlaps with any acgh calls
# (from Dan's cyto tables).

# I've implemented cyto_overlap column. If you see a CNV call that is "concordant" but is FALSE for cyto_overlap
# it means that it matches Dan's manual curated coordinates, but not what the cyto software calls it.

decon_with_acgh <- assign_acgh_status(perrun_decon_calls, acgh_ran_samples)
decon_with_acgh <- assign_non_vscnv(decon_with_acgh, vs_cnv_calls_perrun)
decon_with_acgh <- assign_cyto_overlap(decon_with_acgh, cyto_master_table)

clearcnv_with_acgh <- assign_acgh_status(perrun_clearcnv_calls, acgh_ran_samples)
clearcnv_with_acgh <- assign_non_vscnv(clearcnv_with_acgh, vs_cnv_calls_perrun)
clearcnv_with_acgh <- assign_cyto_overlap(clearcnv_with_acgh, cyto_master_table)


# decon_with_acgh <- decon_with_acgh %>% filter(acgh_status != "tbc")
# clearcnv_with_acgh <- clearcnv_with_acgh %>% filter(acgh_status != "tbc")


decon_with_acgh <- decon_with_acgh %>%
  filter(Quality >= 10)


# ---- Convert calls dataframe to GRanges with Sample info ----
df_to_gr <- function(df) {
  GRanges(
    seqnames = df$Chromosome,
    ranges = IRanges(start = df$Start, end = df$End),
    Sample = df$Sample,
    acgh_status = df$acgh_status
  )
}

decon_gr <- df_to_gr(decon_with_acgh)
clearcnv_gr <- df_to_gr(clearcnv_with_acgh)

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

# ---- Pairwise overlap (Decon vs ClearCNV) ----
decon_clearcnv <- get_overlaps(decon_gr, clearcnv_gr)

# ---- Find tool-specific calls (not shared with the other tool) ----
get_unique_calls <- function(gr1, gr2) {
  hits <- findOverlaps(gr1, gr2)
  hits <- hits[gr1$Sample[queryHits(hits)] == gr2$Sample[subjectHits(hits)]]
  
  unique_idx <- setdiff(seq_along(gr1), unique(queryHits(hits)))
  gr1[unique_idx]
}

decon_only <- get_unique_calls(decon_gr, clearcnv_gr)
decon_only_df <- as.data.frame(decon_only)
clearcnv_only <- get_unique_calls(clearcnv_gr, decon_gr)
clearcnv_only_df <- as.data.frame(clearcnv_only)

# ---- Inspect results ----
decon_clearcnv

sum(decon_only$acgh_status == "concordant")
sum(decon_only$acgh_status == "discordant")
sum(decon_only$acgh_status == "tbc")
sum(decon_only$acgh_status == "inconclusive")

sum(clearcnv_only$acgh_status == "concordant")
sum(clearcnv_only$acgh_status == "discordant")
sum(clearcnv_only$acgh_status == "tbc")
sum(clearcnv_only$acgh_status == "inconclusive")


dim(decon_with_acgh[which(decon_with_acgh$acgh_status == "concordant"), ])[1]
dim(decon_with_acgh[which(decon_with_acgh$acgh_status == "discordant"), ])[1]
dim(decon_with_acgh[decon_with_acgh$acgh_status %in% c("tbc", "inconclusive"), ])[1]

dim(clearcnv_with_acgh[which(clearcnv_with_acgh$acgh_status == "concordant"), ])[1]
dim(clearcnv_with_acgh[which(clearcnv_with_acgh$acgh_status == "discordant"), ])[1]
dim(clearcnv_with_acgh[clearcnv_with_acgh$acgh_status %in% c("tbc", "inconclusive"), ])[1]


# write.csv(decon_with_acgh %>% filter(acgh_status %in% c("tbc", "inconclusive")), "/path/to/project/Aim_2/decon_with_tbc_inconclusive_acgh.csv")
# write.csv(clearcnv_with_acgh %>% filter(acgh_status %in% c("tbc", "inconclusive")), "/path/to/project/Aim_2/clearcnv_with_tbc_inconclusive_acgh.csv")


