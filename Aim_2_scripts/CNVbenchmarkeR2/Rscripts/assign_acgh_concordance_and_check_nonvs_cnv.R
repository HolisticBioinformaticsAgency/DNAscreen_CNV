library(dplyr)
library(GenomicRanges)
library(purrr)
library(tidyr)

#### This script generates the accuracy metrics for all 12 callers
#### and also evaluates each call from each CNV caller and assigns them as 
#### TP (concordant), FP (discordant), or tbc (not in acgh_ran_samples)

## Make sure vscnv_reorder_perrun is created from overlap_calls_from_tools_and_vscnv.R
# validated_res <- read.csv("/path/to/project/Aim_1/CNV_calling/cnvbenchmarker/CNVbenchmarkeR2_test/examples/final_validated_run_val_samples.tsv", sep = "\t")
validated_res <- read.csv("/path/to/project/Aim_1/Validations/acgh_validation_status_files/acgh_round_one_res_for_candidate_cnvs_dan.csv")

## omitting the maybe, likely and "a" because they are CNV calls that Dan cannot be certain if
## concordant or discordant.
# validated_res <- validated_res %>% filter(!Verified_20250730 %in% c("MAYBE", "Likely", "a"))

validated_res$Man..cur..Chr <- gsub("^chr", "", validated_res$Man..cur..Chr)  # remove 'chr' prefix
validated_res$Man..cur..Chr[validated_res$Man..cur..Chr == ""] <- NA            # convert empty strings to NA

validated_res$Var.Seq.Chr <- gsub("^chr", "", validated_res$Var.Seq.Chr)  # remove 'chr' prefix
validated_res$Var.Seq.Chr[validated_res$Var.Seq.Chr == ""] <- NA      

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


# Prepare verified_inconclusive CNVs
verified_inconclusive_cnv <- validated_res %>%
  filter(
    Verified_20250730 %in% c("MAYBE", "Likely", "a")
  ) %>%
  mutate(
    verified_by_dan = Verified_20250730,
    cn = recode(Estimated.CN.of.CNV, "1" = "deletion", "3" = "duplication"),
    
    # Replace missing manual curation coordinates with VarSeq coordinates
    Manual.curation.Start = coalesce(Manual.curation.Start, VarSeq.Start),
    Manual.curation.Stop  = coalesce(Manual.curation.Stop,  VarSeq.Stop),
    Man..cur..Chr         = coalesce(Man..cur..Chr,         Var.Seq.Chr)
  ) %>%
  dplyr::select(
    sample_id_trimmed,
    Gene.Names,
    verified_by_dan,
    Manual.curation.Start,
    Manual.curation.Stop,
    Man..cur..Chr
  )


# This function just assigns acgh_status column for a tool's data and returns enriched data
# discordance is assigned for CNV that does not overlap with any verified==YES CNVs
# Note that for decon, clearcnv and gatk, we checked that there are NO cnv calls 
# (from the 77 array samples) that are not called by vs-cnv that overlaps with any acgh calls 
# (from Dan's cyto tables). So if a CNV call is non-vscnv and is discordant
# it means that the cnv-call is simply not called by acgh as well
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
          val_yes_sub <- verified_yes_cnv %>% filter(sample_id_trimmed == curr_sample)
          val_incon_sub <- verified_inconclusive_cnv %>% filter(sample_id_trimmed == curr_sample)
          
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
              # Helper function to check if CNV overlaps validated CNVs >=50% of ROIs
              check_overlap <- function(val_sub) {
                if (nrow(val_sub) == 0) return(FALSE)
                
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
                  
                  intersect_width <- max(0, min(end(roi_range), end(val_range)) - 
                                           max(start(roi_range), start(val_range)) + 1)
                  frac <- intersect_width / width(roi_range)
                  
                  if (frac >= 0.2) {
                    overlapping_roi_ids <- c(overlapping_roi_ids, mcols(roi_in_call)$roi_id[roi_idx])
                  }
                }
                
                num_overlap <- length(overlapping_roi_ids)
                num_total <- length(roi_in_call)
                
                if (num_total > 0 && (num_overlap / num_total) >= 0.5) TRUE else FALSE
              }
              
              if (check_overlap(val_yes_sub)) {
                "concordant"
              } else if (check_overlap(val_incon_sub)) {
                "inconclusive"
              } else {
                "discordant"
              }
            }
          }
        }
      }
    ) %>%
    ungroup()
}


## The above just assigns discordance to cnv not verified==YES. Function below
## checks if any of the CNV not verified by Dan also appears in the cyto table (called by acgh)
## so they are potential cnvs to look for

gene_whitelist <- c("BRCA1", "BRCA2", "PALB2", "MLH1", "MSH2", "MSH6", "APOB", "PCSK9", "LDLR")

# Prepare cyto table
cyto_master_table <- read.csv(
  "/path/to/project/Aim_1/Validations/stitching_dan_cyto_files/cyto_master_table.tsv", 
  sep = "\t"
) %>%
  mutate(
    sample_id_trimmed = gsub("^.*?_(.*?)_.*$", "\\1", sample_id),
    Chr = sub("^chr", "", Chr)
  ) %>%
  separate_rows(Gene.Names, sep = ",") %>%          # split genes into separate rows
  mutate(Gene.Names = trimws(Gene.Names)) %>%       # trim whitespace
  filter(Gene.Names %in% gene_whitelist)            # keep only whitelist genes


### work on running assign_cyto_overlap on the tool data both validation run and per-run calls

assign_cyto_overlap <- function(tool_data, cyto_master_table) {
  gene_whitelist <- c("BRCA1", "BRCA2", "PALB2", "MLH1", 
                      "MSH2", "MSH6", "APOB", "PCSK9", "LDLR")
  
  # Clean cyto table and derive CNV type
  cyto_clean <- cyto_master_table %>%
    mutate(
      sample_id_trimmed = gsub("^.*?_(.*?)_.*$", "\\1", sample_id),
      Chr = sub("^chr", "", Chr),
      CNV_Type = case_when(
        (Loss != 0 | Deletion != 0) ~ "deletion",
        (Gain != 0 | Amplification != 0) ~ "duplication",
        TRUE ~ NA_character_
      )
    ) %>%
    tidyr::separate_rows(Gene.Names, sep = ",") %>%
    mutate(Gene.Names = trimws(Gene.Names)) %>%
    filter(Gene.Names %in% gene_whitelist, !is.na(CNV_Type))
  
  # Keep all tool calls
  tool_filtered <- tool_data
  
  # Build GRanges with CNV type
  tool_gr <- GenomicRanges::GRanges(
    seqnames = tool_filtered$Chromosome,
    ranges = IRanges::IRanges(start = tool_filtered$Start, end = tool_filtered$End),
    CNV_Type = tool_filtered$CNV_Type,
    sample_id_trimmed = tool_filtered$sample_id_trimmed
  )
  cyto_gr <- GenomicRanges::GRanges(
    seqnames = cyto_clean$Chr,
    ranges = IRanges::IRanges(start = cyto_clean$Start, end = cyto_clean$Stop),
    CNV_Type = cyto_clean$CNV_Type,
    sample_id_trimmed = cyto_clean$sample_id_trimmed
  )
  
  # Find overlaps (without pre-filtering)
  hits <- GenomicRanges::findOverlaps(tool_gr, cyto_gr)
  
  # Default: no overlaps
  tool_filtered$cyto_overlap <- FALSE
  
  if (length(hits) > 0) {
    for (i in seq_along(hits)) {
      q <- queryHits(hits)[i]
      s <- subjectHits(hits)[i]
      
      # Require CNV type *and* sample to match
      if (tool_gr$CNV_Type[q] == cyto_gr$CNV_Type[s] &&
          tool_gr$sample_id_trimmed[q] == cyto_gr$sample_id_trimmed[s]) {
        tool_filtered$cyto_overlap[q] <- TRUE
      }
    }
  }
  
  return(tool_filtered)
}



# Function to assign 'non_vscnv' status
assign_non_vscnv <- function(tool_data, vs_cnv_calls_perrun) {
  tool_data %>%
    mutate(sample_id_trimmed = sub(".*-(DNS\\d+)_S\\d+.*", "\\1", Sample)) %>%
    rowwise() %>%
    mutate(
      non_vscnv = {
        # CNV call from the tool
        call_gr <- GRanges(seqnames = Chromosome, ranges = IRanges(start = Start, end = End))
        curr_sample <- sample_id_trimmed
        curr_gene <- Gene
        
        # VS-CNV calls for this sample
        vs_sub <- vs_cnv_calls_perrun %>%
          filter(sample_id_trimmed == curr_sample)
        
        if (nrow(vs_sub) == 0) {
          # No VS-CNV calls for this sample → call is unique and non-vscnv
          "yes"
        } else {
          # ROIs for this gene
          rois <- bed_table_clean %>%
            filter(gene == curr_gene)
          
          if (nrow(rois) == 0) {
            "yes" # no ROIs in the nine panel gene → consider unique and non-vscnv
          } else {
            roi_gr <- GRanges(
              seqnames = rois$chrom,
              ranges = IRanges(start = rois$start, end = rois$stop)
            )
            
            # Find ROIs overlapping with the tool's CNV call
            hits <- findOverlaps(roi_gr, call_gr)
            roi_in_call <- roi_gr[queryHits(hits)]
            
            if (length(roi_in_call) == 0) {
              "yes" # no ROIs overlapped (intronic?) → consider unique and non-vscnv
            } else {
              # For each VS-CNV call, find its overlapping ROIs
              vs_rois_list <- list()
              
              for (i in 1:nrow(vs_sub)) {
                vs_call <- vs_sub[i, ]
                vs_gr <- GRanges(
                  seqnames = vs_call$VS_Call_chr,
                  ranges = IRanges(start = vs_call$VS_Call_start, end = vs_call$VS_Call_end)
                )
                
                # Find ROIs overlapping with this VS-CNV call
                vs_hits <- findOverlaps(roi_gr, vs_gr)
                vs_roi_in_call <- roi_gr[queryHits(vs_hits)]
                
                # Store the genomic ranges for this VS-CNV call
                vs_rois_list[[i]] <- vs_roi_in_call
              }
              
              # Combine all ROIs from all VS-CNV calls
              all_vs_rois <- unique(do.call("c", vs_rois_list))
              
              if (length(all_vs_rois) == 0) {
                "yes" # no ROIs in VS-CNV calls → unique
              } else {
                # Count how many ROIs in the tool's call overlap with VS-CNV ROIs
                overlap_count <- 0
                total_rois <- length(roi_in_call)
                
                for (j in 1:total_rois) {
                  tool_roi <- roi_in_call[j]
                  
                  # Check if this ROI overlaps with any VS-CNV ROI
                  overlaps <- findOverlaps(tool_roi, all_vs_rois)
                  
                  if (length(overlaps) > 0) {
                    # For each overlap, check if it meets the 50% threshold
                    for (k in 1:length(overlaps)) {
                      vs_roi <- all_vs_rois[subjectHits(overlaps)[k]]
                      
                      # Calculate overlap percentage
                      tool_range <- ranges(tool_roi)
                      vs_range <- ranges(vs_roi)
                      intersect_width <- max(0, min(end(tool_range), end(vs_range)) - 
                                               max(start(tool_range), start(vs_range)) + 1)
                      overlap_frac <- intersect_width / width(tool_range)
                      
                      if (overlap_frac >= 0.5) {
                        overlap_count <- overlap_count + 1
                        break # count this ROI only once even if it overlaps multiple VS-CNV ROIs
                      }
                    }
                  }
                }
                
                # Check if more than half of the tool's ROIs overlap with VS-CNV ROIs
                if (overlap_count / total_rois > 0) {
                  "no" # more than half of ROIs match VS-CNV → not unique
                } else {
                  "yes" # unique call
                }
              }
            }
          }
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


harmonised_list_with_nonvscnv <- imap(harmonised_list_with_acgh, ~ {
  enriched <- assign_non_vscnv(.x$data, vs_cnv_calls_perrun)
  list(data = enriched, quality_col = .x$quality_col)
})

harmonised_list_with_nonvscnv <- imap(harmonised_list_with_nonvscnv, ~ {
  enriched <- assign_cyto_overlap(.x$data, cyto_master_table)
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

## Summary of number of calls that are vscnv and non_vscnv

callers_to_use <- setdiff(
  names(harmonised_list_with_nonvscnv),
  c("VS-CNV_perrun", "VS-CNV_val")
)

summary_table <- lapply(callers_to_use, function(caller) {
  
  caller_data <- harmonised_list_with_nonvscnv[[caller]]$data
  caller_data <- caller_data %>% filter(acgh_status != "tbc")
  
  x <- caller_data$non_vscnv
  
  total_calls <- sum(!is.na(x))
  non_vscnv_calls <- sum(x == "yes", na.rm = TRUE)
  percent_non_vscnv <- (non_vscnv_calls / total_calls) * 100
  
  data.frame(
    caller = caller,
    total_calls = total_calls,
    non_vscnv_calls = non_vscnv_calls,
    percent_non_vscnv = percent_non_vscnv
  )
})

summary_table <- do.call(rbind, summary_table)

overall_total_calls <- sum(summary_table$total_calls)
overall_non_vscnv_calls <- sum(summary_table$non_vscnv_calls)
overall_percent_non_vscnv <- (overall_non_vscnv_calls / overall_total_calls) * 100

summary_table <- rbind(
  summary_table,
  data.frame(
    caller = "ALL_CALLERS",
    total_calls = overall_total_calls,
    non_vscnv_calls = overall_non_vscnv_calls,
    percent_non_vscnv = overall_percent_non_vscnv
  )
)



# get all_metrics from the new enriched list
# all_metrics <- imap_dfr(harmonised_list_with_acgh, ~ evaluate_tool(list(data = .x$data), .y, acgh_ran_samples))

# write.table(all_metrics, "output/summary/results_val_round_one.tsv",sep = "\t", row.names = FALSE)



