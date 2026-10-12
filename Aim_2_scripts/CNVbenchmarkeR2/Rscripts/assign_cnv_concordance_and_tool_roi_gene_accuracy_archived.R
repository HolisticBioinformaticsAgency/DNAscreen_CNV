library(dplyr)
library(GenomicRanges)
library(purrr)
library(stringr)

#### This script generates the accuracy metrics for all 12 callers
#### and also evaluates each call from each CNV caller and assigns them as 
#### TP (concordant), FP (discordant), or tbc (not in acgh_ran_samples)

## Make sure vscnv_reorder_perrun is created from overlap_calls_from_tools_and_vscnv.R
# validated_res <- read.csv("/path/to/project/Aim_1/CNV_calling/cnvbenchmarker/CNVbenchmarkeR2_test/examples/final_validated_run_val_samples.tsv", sep = "\t")
validated_res <- read.csv("/path/to/project/Aim_1/Validations/acgh_validation_status_files/acgh_round_one_res_for_candidate_cnvs_dan.csv")

## omitting the maybe, likely and "a" because they are CNV calls that Dan cannot be certain if
## concordant or discordant
validated_res <- validated_res %>% filter(!Verified_20250730 %in% c("MAYBE", "Likely", "a"))

validated_res$Man..cur..Chr <- gsub("^chr", "", validated_res$Man..cur..Chr)  # remove 'chr' prefix
validated_res$Man..cur..Chr[validated_res$Man..cur..Chr == ""] <- NA            # convert empty strings to NA

validated_res[] <- lapply(validated_res, function(x) {
  if (is.character(x)) x[x == ""] <- NA
  return(x)
})


# bed_table <- read.csv("/path/to/project/Aim_1/CNV_calling/exome_depth/bed_files/9genes_25bp.fix.sorted.bed", 
#                       sep = "\t", col.names = c("chrom", "start", "stop", "roi", "idk", "strand"))
bed_table <- read.csv("/path/to/project/Aim_2/cnvbenchmarker/CNVbenchmarkeR2_test/bed_file/3427651_Covered_v3.bed", 
                      sep = "\t", col.names = c("chrom", "start", "stop", "roi"))

# Clean BED table
bed_table_clean <- bed_table %>%
  mutate(chrom = gsub("^chr", "", chrom),
         gene = sub("_cds_.*", "", roi))

acgh_ran_samples <- unique(validated_res$sample_id_trimmed)

# Prepare verified_yes CNVs
verified_yes_cnv <- validated_res %>%
  filter(
    Verified_20250730 == "YES"
  ) %>%
 mutate(
    verified_by_dan = Verified_20250730,
    cn = recode(Estimated.CN.of.CNV, "1" = "deletion", "3" = "duplication"),
  ) %>%
  dplyr::select(sample_id_trimmed, Gene.Names, verified_by_dan, Manual.curation.Start, Manual.curation.Stop, Man..cur..Chr)

# Convert manual curation coordinates to GRanges
verified_yes_cnv_gr <- GRanges(
  seqnames = verified_yes_cnv$Man..cur..Chr,
  ranges = IRanges(start = verified_yes_cnv$Manual.curation.Start, end = verified_yes_cnv$Manual.curation.Stop)
)

# Perform overlaps in order to assign gene names to verified_yes CNV (for gene-level accuracy calculation)
# Read the gene coordinates
gene_data <- read.table("/path/to/project/Aim_1/CNV_calling/CI_plots/dnascreen_genes.txt", header = FALSE, col.names = c("chr", "start", "end", "Gene.Names"))
gene_data$chr <- str_replace(gene_data$chr, "^chr", "")
# Convert gene_data to GRanges
gene_ranges <- GRanges(
  seqnames = gene_data$chr,
  ranges = IRanges(start = gene_data$start, end = gene_data$end),
  gene_name = gene_data$Gene.Names
)
# gene_ranges from overlap_calls_from_tools_and_vscnv.R
overlaps <- findOverlaps(verified_yes_cnv_gr, gene_ranges)

# Extract overlapping gene names
verified_yes_cnv$Gene.from.coords <- NA  # placeholder
verified_yes_cnv$Gene.from.coords[queryHits(overlaps)] <- gene_ranges$gene_name[subjectHits(overlaps)]
verified_yes_cnv$Gene <- verified_yes_cnv$Gene.from.coords

# Generate all sample-gene combinations
all_sample_gene_combos <- expand.grid(
  sample_id_trimmed = unique(verified_yes_cnv$sample_id_trimmed),
  gene = unique(bed_table_clean$gene),
  stringsAsFactors = FALSE
)


### Function to calculate accuracy metrics for VS-CNV like how CNVBenchmarkeR2 does it! ###

calc_metrics <- function(tp, fp, fn, tn, algorithm = "vs_cnv", level) {
  total <- tp + fp + fn + tn
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA
  ppv         <- if ((tp + fp) > 0) tp / (tp + fp) else NA  # precision
  npv         <- if ((tn + fn) > 0) tn / (tn + fn) else NA
  f1_score    <- if ((ppv + sensitivity) > 0) 2 * ppv * sensitivity / (ppv + sensitivity) else NA
  
  mcc_num <- (tp * tn - fp * fn)
  mcc_den <- sqrt(as.numeric(tp + fp) * as.numeric(tp + fn) * 
                    as.numeric(tn + fp) * as.numeric(tn + fn))
  mcc <- if (!is.na(mcc_den) && mcc_den > 0) mcc_num / mcc_den else NA
  
  po <- (tp + tn) / total
  pe <- (((tp + fp) * (tp + fn)) + ((tn + fn) * (tn + fp))) / (total^2)
  kappa <- if ((1 - pe) != 0) (po - pe) / (1 - pe) else NA
  
  accuracy <- (tp + tn) / total
  
  tibble::tibble(
    TP = tp,
    TN = tn,
    FP = fp,
    FN = fn,
    total = total,
    sensitivity = round(sensitivity, 4),
    specificity = round(specificity, 4),
    PPV = round(ppv, 4),
    NPV = round(npv, 4),
    F1 = round(f1_score, 4),
    MCC = round(mcc, 4),
    kappa = round(kappa, 4),
    accuracy = round(accuracy, 4),
    algorithm = algorithm,
    level = level
  )
}

# This function just assigns acgh_status column for a tool's data and returns enriched data
assign_acgh_status <- function(tool_data, acgh_ran_samples) {
  tool_data %>%
    mutate(sample_id_trimmed = sub(".*-(DNS\\d+)_S\\d+.*", "\\1", Sample)) %>%
    rowwise() %>%
    mutate(
      acgh_status = {
        if (!(sample_id_trimmed %in% acgh_ran_samples)) {
          "tbc"
        } else {
          # CNV call from the tool
          call_gr <- GRanges(seqnames = Chromosome,
                             ranges = IRanges(start = Start, end = End))
          curr_sample <- sample_id_trimmed
          curr_gene <- Gene
          
          # aCGH validated CNVs for this sample
          val_sub <- verified_yes_cnv %>%
            filter(sample_id_trimmed == curr_sample)
          
          if (nrow(val_sub) == 0) {
            "discordant"
          } else {
            # ROIs for this gene
            rois <- bed_table_clean %>%
              filter(gene == curr_gene) %>%
              mutate(roi_id = roi)
            
            if (nrow(rois) == 0) {
              "discordant"
            } else {
              roi_gr <- GRanges(
                seqnames = rois$chrom,
                ranges = IRanges(start = rois$start, end = rois$stop),
                roi_id = rois$roi_id
              )
              
              # Keep only ROIs overlapping with CNV call
              hits <- findOverlaps(roi_gr, call_gr)
              roi_in_call <- roi_gr[queryHits(hits)]
              
              if (length(roi_in_call) == 0) {
                "discordant"
              } else {
                # aCGH GRanges for this sample
                val_gr <- GRanges(
                  seqnames = val_sub$Man..cur..Chr,
                  ranges = IRanges(start = val_sub$Manual.curation.Start,
                                   end = val_sub$Manual.curation.Stop)
                )
                
                overlaps_df <- as.data.frame(findOverlaps(roi_in_call, val_gr))
                overlapping_roi_ids <- c()
                
                for (i in seq_len(nrow(overlaps_df))) {
                  roi_idx <- overlaps_df$queryHits[i]
                  val_idx <- overlaps_df$subjectHits[i]
                  
                  roi_range <- ranges(roi_in_call)[roi_idx]
                  val_range <- ranges(val_gr)[val_idx]
                  
                  intersect_width <- max(0, min(end(roi_range), end(val_range)) - max(start(roi_range), start(val_range)) + 1)
                  frac <- intersect_width / width(roi_range)
                  
                  if (frac >= 0.2) {
                    overlapping_roi_ids <- c(overlapping_roi_ids, mcols(roi_in_call)$roi_id[roi_idx])
                  }
                }

                num_overlap <- length((overlapping_roi_ids))
                num_total <- length(roi_in_call)

                if (num_total > 0 && (num_overlap / num_total) >= 0.5) {
                  "concordant"
                } else {
                  "discordant"
                }
              }
            }
          }
        }
      }
    ) %>%
    ungroup()
}


# Function to assign 'non_vscnv' status
assign_non_vscnv <- function(tool_data, vs_cnv_calls_perrun) {
  # Prepare VS-CNV GRanges
  vs_gr <- GRanges(
    seqnames = vs_cnv_calls_perrun$VS_Call_chr,
    ranges = IRanges(start = vs_cnv_calls_perrun$VS_Call_start,
                     end = vs_cnv_calls_perrun$VS_Call_end),
    sample = vs_cnv_calls_perrun$Sample
  )
  
  tool_data %>%
    mutate(sample_id_trimmed = sub(".*-(DNS\\d+)_S\\d+.*", "\\1", Sample)) %>%
    rowwise() %>%
    mutate(
      non_vscnv = {
        call_gr <- GRanges(
          seqnames = Chromosome,
          ranges = IRanges(start = Start, end = End),
          sample = sample_id_trimmed
        )
        
        # Find overlaps with VS-CNV for the same sample
        hits <- findOverlaps(call_gr, vs_gr)
        hits_df <- as.data.frame(hits)
        
        # Filter by same sample
        hits_df <- hits_df[mcols(call_gr)$sample[hits_df$queryHits] == 
                             mcols(vs_gr)$sample[hits_df$subjectHits], ]
        
        if (nrow(hits_df) == 0) {
          "yes"  # no overlapping VS-CNV → unique call
        } else {
          # check fraction overlap > 50%
          overlaps_frac <- sapply(seq_len(nrow(hits_df)), function(i) {
            q_range <- ranges(call_gr)[hits_df$queryHits[i]]
            s_range <- ranges(vs_gr)[hits_df$subjectHits[i]]
            intersect_width <- max(0, min(end(q_range), end(s_range)) - max(start(q_range), start(s_range)) + 1)
            intersect_width / width(q_range)
          })
          
          if (any(overlaps_frac >= 0.5)) "no" else "yes"
        }
      }
    ) %>%
    ungroup()
}


# Create harmonised_list_with_acgh: add acgh_status to each tool's data
harmonised_list_with_acgh <- imap(harmonised_list, ~ {
  enriched <- assign_acgh_status(.x$data, acgh_ran_samples)
  # Return a modified tool entry with enriched data and keep other info
  list(data = enriched, quality_col = .x$quality_col)
})



# Main loop function for one tool
evaluate_tool <- function(tool_entry, tool_name, acgh_ran_samples) {
  cat("Evaluating:", as.character(tool_name), "\n")
  
  tool_data <- tool_entry$data
  
  # Trim sample ID and select relevant columns
  tool_calls <- tool_data %>%
    mutate(sample_id_trimmed = sub(".*-(DNS\\d+)_S\\d+.*", "\\1", Sample)) %>%
    dplyr::select(sample_id_trimmed, Gene, Chromosome, Start, End)
  
  # Assign 'acgh_status' for each CNV call in tool_calls
  tool_calls <- assign_acgh_status(tool_entry$data, acgh_ran_samples) %>%
    dplyr::select(sample_id_trimmed, Gene, Chromosome, Start, End, acgh_status)
  
  # For ROI-level: sample × individual ROI
  roi_combos <- expand.grid(
    sample_id_trimmed = acgh_ran_samples,
    roi_index = seq_len(nrow(bed_table_clean)),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      gene = bed_table_clean$gene[roi_index],
      chrom = bed_table_clean$chrom[roi_index],
      start = bed_table_clean$start[roi_index],
      stop = bed_table_clean$stop[roi_index],
      roi_id = bed_table_clean$roi[roi_index]
    )
  
  
  ## ROI-level metrics
  tp_roi <- 0; fp_roi <- 0; fn_roi <- 0; tn_roi <- 0
  
  ### only count as FN if overlap <20%...
  
  roi_combos <- expand.grid(
    sample_id_trimmed = acgh_ran_samples,
    roi_index = seq_len(nrow(bed_table_clean)),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      gene = bed_table_clean$gene[roi_index],
      chrom = bed_table_clean$chrom[roi_index],
      start = bed_table_clean$start[roi_index],
      stop = bed_table_clean$stop[roi_index],
      roi_id = bed_table_clean$roi[roi_index]
    ) %>%
    filter(gene != "ATM")  # omit ATM
  
  # Create one GRanges object for all ROI combos
  roi_gr <- GRanges(
    seqnames = roi_combos$chrom,
    ranges = IRanges(start = roi_combos$start, end = roi_combos$stop),
    sample = roi_combos$sample_id_trimmed,
    gene = roi_combos$gene,
    roi_id = roi_combos$roi_id
  )
  
  # Create GRanges for validated CNVs
  validated_gr <- GRanges(
    seqnames = verified_yes_cnv$Man..cur..Chr,
    ranges = IRanges(start = verified_yes_cnv$Manual.curation.Start,
                     end = verified_yes_cnv$Manual.curation.Stop),
    sample = verified_yes_cnv$sample_id_trimmed,
    gene = verified_yes_cnv$Gene
  )
  
  # Create GRanges for tool calls (exclude "tbc")
  tool_calls_filtered <- tool_calls %>% filter(acgh_status != "tbc")
  call_gr <- GRanges(
    seqnames = tool_calls_filtered$Chromosome,
    ranges = IRanges(start = tool_calls_filtered$Start, end = tool_calls_filtered$End),
    sample = tool_calls_filtered$sample_id_trimmed,
    gene = tool_calls_filtered$Gene
  )
  
  # Function to flag ROI overlaps (>=20%) with another GRanges
  get_roi_hits <- function(roi_gr, subject_gr) {
    olaps <- findOverlaps(roi_gr, subject_gr)
    if (length(olaps) == 0) return(logical(length(roi_gr)))
    
    # Only keep overlaps where samples and genes match
    olaps_df <- as.data.frame(olaps)
    match_sample <- mcols(roi_gr)$sample[olaps_df$queryHits] == mcols(subject_gr)$sample[olaps_df$subjectHits]
    match_gene <- mcols(roi_gr)$gene[olaps_df$queryHits] == mcols(subject_gr)$gene[olaps_df$subjectHits]
    olaps_df <- olaps_df[match_sample & match_gene, ]
    
    if (nrow(olaps_df) == 0) return(logical(length(roi_gr)))
    
    # Compute intersection fraction
    query_ranges <- ranges(roi_gr)[olaps_df$queryHits]
    subject_ranges <- ranges(subject_gr)[olaps_df$subjectHits]
    intersection_width <- pmax(0, pmin(end(query_ranges), end(subject_ranges)) - pmax(start(query_ranges), start(subject_ranges)) + 1)
    roi_width <- width(query_ranges)
    frac <- intersection_width / roi_width
    
    hit_flags <- logical(length(roi_gr))
    hit_flags[olaps_df$queryHits[frac >= 0.2]] <- TRUE ## Note that the overlap is kinda redundant for call_hits as most tools will include the entire ROI as part of the call
    return(hit_flags)
  }
  
  # ROI-level overlap calls
  val_hits <- get_roi_hits(roi_gr, validated_gr)
  call_hits <- get_roi_hits(roi_gr, call_gr)
  
  # Confusion matrix
  tp_roi <- sum(val_hits & call_hits)
  fp_roi <- sum(!val_hits & call_hits)
  fn_roi <- sum(val_hits & !call_hits)
  tn_roi <- sum(!val_hits & !call_hits)
  ### Note that in this logic, any ROI that <20% overlap with the aCGH call will not be assessed as we do not expect RD tools
  ### to reliably detect it.
  
  roi_metrics <- calc_metrics(tp_roi, fp_roi, fn_roi, tn_roi, algorithm = tool_name, level = "roi")
  
  ## Gene-level metrics
  tp_gene <- 0; fp_gene <- 0; fn_gene <- 0; tn_gene <- 0
  
  # Get unique gene/sample combos to evaluate
  gene_combos <- expand.grid(
    sample_id_trimmed = acgh_ran_samples,
    gene = setdiff(unique(bed_table_clean$gene), "ATM"),
    stringsAsFactors = FALSE
  )
  
  # Flag presence of validated CNV per sample-gene
  val_gene_status <- verified_yes_cnv %>%
    filter(sample_id_trimmed %in% acgh_ran_samples) %>%
    distinct(sample_id_trimmed, gene = Gene) %>%
    mutate(val_has_cnv = TRUE)
  
  # Flag presence of tool call per sample-gene (excluding tbc)
  call_gene_status <- tool_calls %>%
    filter(acgh_status != "tbc") %>%
    distinct(sample_id_trimmed, gene = Gene) %>%
    mutate(call_has_cnv = TRUE)
  
  # Combine with all gene/sample combos
  gene_status <- gene_combos %>%
    left_join(val_gene_status, by = c("sample_id_trimmed", "gene")) %>%
    left_join(call_gene_status, by = c("sample_id_trimmed", "gene")) %>%
    mutate(
      val_has_cnv = ifelse(is.na(val_has_cnv), FALSE, val_has_cnv),
      call_has_cnv = ifelse(is.na(call_has_cnv), FALSE, call_has_cnv)
    )
  
  # Compute confusion matrix
  tp_gene <- sum(gene_status$val_has_cnv & gene_status$call_has_cnv)
  fp_gene <- sum(!gene_status$val_has_cnv & gene_status$call_has_cnv)
  fn_gene <- sum(gene_status$val_has_cnv & !gene_status$call_has_cnv)
  tn_gene <- sum(!gene_status$val_has_cnv & !gene_status$call_has_cnv)
  
  gene_metrics <- calc_metrics(tp_gene, fp_gene, fn_gene, tn_gene,
                               algorithm = tool_name, level = "gene")
  
  result <- bind_rows(roi_metrics, gene_metrics)
  cat("Tool:", tool_name, "- Total ROIs evaluated:", tp_roi + fp_roi + fn_roi + tn_roi, "\n")
  cat("Tool:", tool_name, "- Total Genes evaluated:", tp_gene + fp_gene + fn_gene + tn_gene, "\n")
  return(result)
}

# New function to return both metrics and enriched data
evaluate_tool_and_return_data <- function(tool_entry, tool_name, acgh_ran_samples) {
  enriched_data <- assign_acgh_status(tool_entry$data, acgh_ran_samples)
  metrics <- evaluate_tool(list(data = enriched_data), tool_name, acgh_ran_samples)
  list(metrics = metrics, data = enriched_data)
}



# get all_metrics from the new enriched list
all_metrics <- imap_dfr(harmonised_list_with_acgh, ~ evaluate_tool(list(data = .x$data), .y, acgh_ran_samples))

# write.table(all_metrics, "output/summary/results_val_round_one.tsv",sep = "\t", row.names = FALSE)



