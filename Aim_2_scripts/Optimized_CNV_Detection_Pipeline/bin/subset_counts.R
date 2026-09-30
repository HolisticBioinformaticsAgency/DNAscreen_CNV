#!/usr/bin/env Rscript
# Subset a DECoN ReadInBams.R RData to the retained (non-outlier) BAMs so
# DECoN can reuse the stage-1 counts instead of re-reading every BAM.
# Usage: subset_counts.R <in.RData> <retained_bams.txt> <out.RData>

args <- commandArgs(trailingOnly = TRUE)
load(args[1])                                   # counts, bams, bed.file, sample.names, fasta
keep <- readLines(args[2])
keep <- keep[nzchar(keep)]
idx  <- match(keep, bams)
if (anyNA(idx)) stop("Retained BAMs not found in counts RData: ", paste(keep[is.na(idx)], collapse = ", "))

n_before     <- length(bams)
counts       <- counts[, c(1:5, 5 + idx)]
bams         <- bams[idx]
sample.names <- sample.names[idx]
save(counts, bams, bed.file, sample.names, fasta, file = args[3])
message(sprintf("Kept %d of %d BAMs", length(idx), n_before))
