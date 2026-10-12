library(ExomeDepth)
library(tools)
library(dplyr)
library(ggplot2)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(tidyr)
library(biomaRt)
library(patchwork)
library(grid)
library(stringr)
library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)
library(ggrepel)


# ---------- Parameters ----------
directory_bam <- "/path/to/project/Aim_1/Validations/bam_and_vcf"
directory_vcf <- "/path/to/project/Aim_1/Validations/bam_and_vcf"

### loading acgh validation data to separate concordant and discordant calls


## Make sure vscnv_reorder_perrun is created from overlap_calls_from_tools_and_vscnv.R
# validated_res <- read.csv("/path/to/project/Aim_1/CNV_calling/cnvbenchmarker/CNVbenchmarkeR2_test/examples/final_validated_run_val_samples.tsv", sep = "\t")
validated_res <- read.csv("/path/to/project/Aim_1/Validations/acgh_validation_status_files/acgh_round_one_confirmed_cnvs_dan.csv")

validated_res$Man..cur..Chr <- gsub("^chr", "", validated_res$Man..cur..Chr)  # remove 'chr' prefix
validated_res$Man..cur..Chr[validated_res$Man..cur..Chr == ""] <- NA            # convert empty strings to NA

validated_res[] <- lapply(validated_res, function(x) {
  if (is.character(x)) x[x == ""] <- NA
  return(x)
})

validated_res_clean <- validated_res %>%
  rename(
    Sample = sample,
    Chromosome = VS_Call_chr,
    Start = VS_Call_start,
    End = VS_Call_end,
    Gene = Gene.Names,
    CNV_Type = Type.of.CNV
  ) %>%
  arrange(adjusted_p_value) %>% 
  mutate(Sample = sub("\\.hq\\.sorted\\.marked.*", "", Sample),
         Chromosome = gsub("^chr", "", Chromosome),
         CNV_Type = recode(CNV_Type,
                           "Duplicate" = "duplication",
                           "Het Deletion" = "deletion")) %>%
  filter(!is.na(Gene)) %>% ## Filter out any CNVs in ATM ROI
  dplyr::select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

# ---------- Rest of the script ----------

# File paths
read_count_dir <- "/path/to/project/Aim_1/CNV_calling/exome_depth/RData"

target_sample_dots <- gsub("-", ".", target_sample)

ref_fa = '/path/to/project/Aim_1/CNV_calling/exome_depth/Homo_sapiens_assembly38.fasta'
bed_name = "9genes_25bp.fix.sorted"
# bed_name = "3427651_Covered"
bed_file = paste0('/path/to/project/Aim_1/CNV_calling/exome_depth/bed_files/', 
                  bed_name, ".bed")



plot_cnv_for_sample <- function(cnv_row, read_count_dir, bed_name, bed_file, 
                                plot_types = c("read_ratio_plot", "exon_plot", "raw_depth_plot"),
                                gc_correction = TRUE, save_folder) {
  # Extract info from CNV row
  target_sample <- paste0(cnv_row$Sample, ".hq.sorted.marked.bam")
  target_sample_dots <- gsub("-", ".", target_sample)
  user_chr <- as.numeric(cnv_row$Chromosome)
  user_start <- as.numeric(cnv_row$Start)
  user_end <- as.numeric(cnv_row$End)
  cnv_start <- as.numeric(cnv_row$Start)
  cnv_end <- as.numeric(cnv_row$End)
  gene_symbol <- cnv_row$Gene
  cnv_type <- cnv_row$CNV_Type
  dnascreen_run <- cnv_row$run  # Make sure this column is in your df
  
  if (is.na(gene_symbol)) {
    message("Skipping: gene_symbol is NA. Likely a CNV called in ATM which is covering one small variant only")
    return(NULL)
  }
  
  # Coordinates dictionary (gene to coordinates)
  gene_coords <- list(
    BRCA1 = c(17, 43040241, 43129417),
    BRCA2 = c(13, 32311269, 32404506),
    PALB2 = c(16, 23601257, 23643217),
    MLH1 = c(3, 36990475, 37053720),
    MSH2 = c(2, 47399058, 47487236),
    MSH6 = c(2, 47781888, 47808146),
    APOB = c(2, 20999296, 21046205),
    LDLR = c(19, 11087245, 11136037),
    PCSK9 = c(1, 55038282, 55066118)
  )
  
  if (gene_symbol %in% names(gene_coords)) {
    coords <- gene_coords[[gene_symbol]]
    user_chr <- coords[1]
    user_start <- coords[2]
    user_end <- coords[3]
  } else {
    warning(paste("Gene symbol", gene_symbol, "not in coordinate dictionary"))
  }
  
  # You would place all the plotting code from your big script here, wrapped with `if ("plot_type" %in% plot_types)`
  if ("read_ratio_plot" %in% plot_types) {
    print(paste0("plotting read ratio plot for ", target_sample, "_", gene_symbol))
    
    read_count_rdata  <- paste0(read_count_dir, "/run", dnascreen_run, "_bams_", bed_name, ".RData")
    
    if (file.exists(read_count_rdata)){
      load(read_count_rdata)
    } else {
      print("RData file not found. Terminating.")
      stop(paste("Missing file:", read_count_rdata))
      # my.count <- getBamCounts(bed.file = bed_file,
      #                          bam.files = bam_files,
      #                          include.chr = FALSE,
      #                          referenceFasta = ref_fa)
      # save(my.count, file = paste0("/path/to/project/Aim_1/CNV_calling/exome_depth/RData/run", dnascreen_run, "_bams_", bed_name ,".RData"))
    }
    
    
    my.count.dafr <- as(my.count, 'data.frame')
    
    # Get sample ID columns (assuming samples start from column 6)
    sample_ids <- colnames(my.count.dafr)[6:ncol(my.count.dafr)]
    
    # Extract count matrix and ensure numeric type
    count_matrix <- my.count.dafr[, sample_ids]
    count_matrix <- as.data.frame(lapply(count_matrix, as.numeric))
    
    ### Function to apply gc loess normalization
    gc_loess_normalize_wide <- function(read_count_df, gc_content_vector) {
      
      # Ensure GC content is same length as rows
      stopifnot(length(gc_content_vector) == nrow(read_count_df))
      
      # Add GC content column
      read_count_df <- read_count_df %>%
        mutate(region = row_number(), gc_content = gc_content_vector)
      
      # Pivot to long format
      long_df <- read_count_df %>%
        pivot_longer(
          cols = -c(region, gc_content),
          names_to = "sample",
          values_to = "read_count"
        )
      
      # Apply LOESS normalization per sample
      normalized_long_df <- long_df %>%
        group_by(sample) %>%
        mutate(
          loess_fit = loess(read_count ~ gc_content, span = 0.75)$fitted,
          read_count_gc_normalized = read_count / loess_fit * median(read_count, na.rm = TRUE)
        ) %>%
        ungroup()
      
      # Return to wide format
      normalized_wide <- normalized_long_df %>%
        dplyr::select(region, sample, read_count_gc_normalized) %>%
        pivot_wider(names_from = sample, values_from = read_count_gc_normalized) %>%
        arrange(region) %>%
        dplyr::select(-region)
      
      return(normalized_wide)
    }
    
    if (gc_correction == "TRUE") {
      count_matrix <- gc_loess_normalize_wide(count_matrix, my.count$GC)
    }
    
    # Normalize each row by total sample read count (in counts per thousand)
    # Calculate total reads per sample (column-wise sum)
    sample_totals <- colSums(count_matrix, na.rm = TRUE)
    
    # Normalize each sample by its total read count (Counts Per Thousand)
    normalized_counts <- sweep(count_matrix, 2, sample_totals, FUN = "/") * 1e3
    
    normalized_df <- as.data.frame(normalized_counts)
    colnames(normalized_df) <- sample_ids
    
    my.count.dafr[, sample_ids] <- normalized_df
    
    my.count.dafr$chromosome <- gsub("^chr", "", my.count.dafr$chromosome)
    my.count.dafr$chromosome <- as.numeric(my.count.dafr$chromosome)
    
    
    # Filter regions by user-specified coordinates
    filtered_counts <- my.count.dafr %>%
      dplyr::filter(chromosome == user_chr & start >= user_start & end <= user_end)
    
    if (nrow(filtered_counts) == 0) {
      stop("No regions found in the specified coordinates.")
    }
    
    # Extract sample columns and ensure numeric
    sample_ids <- colnames(filtered_counts)[6:ncol(filtered_counts)]
    count_matrix <- filtered_counts[, sample_ids] %>%
      mutate_all(as.numeric)
    
    # # --- Compute Average Percent Difference for each sample ---
    
    # Function to calculate percent difference
    calculate_avg_percent_diff <- function(sample_col, reference_matrix) {
      reference_median <- apply(reference_matrix, 1, median, na.rm = TRUE)
      percent_diff <- abs(sample_col - reference_median) / reference_median * 100
      mean(percent_diff, na.rm = TRUE)
    }
    
    # Create a matrix of just numeric sample data (already normalized)
    normalized_matrix <- my.count.dafr[, sample_ids]
    
    # Initialize APD vector
    avg_percent_diff_vec <- c()
    
    for (sample in sample_ids) {
      reference_matrix <- normalized_matrix[, setdiff(sample_ids, sample)]
      avg_percent_diff_vec[sample] <- calculate_avg_percent_diff(normalized_matrix[[sample]], reference_matrix)
    }
    
    # --- Select Top 10 reference samples with the lowest APD ---
    
    # Sort APD values (excluding target_sample_dots)
    avg_percent_diff_vec_filtered <- avg_percent_diff_vec[names(avg_percent_diff_vec) != target_sample_dots]
    top20_samples <- names(sort(avg_percent_diff_vec_filtered))[1:20]
    
    # Extract region-level read counts for target sample
    filtered_counts <- my.count.dafr %>%
      filter(chromosome == user_chr & start >= user_start & end <= user_end)
    
    if (nrow(filtered_counts) == 0) {
      stop("No regions found in the specified coordinates.")
    }
    
    # Extract sample columns and ensure numeric
    sample_ids <- colnames(filtered_counts)[6:ncol(filtered_counts)]
    count_matrix <- filtered_counts[, sample_ids] %>%
      mutate_all(as.numeric)
    
    # Target sample vector
    target_vector <- count_matrix[[target_sample_dots]]
    
    # Region medians from top 10 reference samples
    region_medians_top20 <- apply(count_matrix[, top20_samples], 1, median, na.rm = TRUE)
    
    # Read ratio
    read_ratio_vector <- target_vector / region_medians_top20
    
    # Prepare data for plotting
    plot_df <- filtered_counts %>%
      mutate(read_ratio = read_ratio_vector) %>%
      dplyr::select(start, end, read_ratio, exon) %>%
      dplyr::arrange(start)
    
    ## For the continuous line
    # Initialize stepwise plot data
    step_df <- data.frame()
    
    for (i in 1:(nrow(plot_df) - 1)) {
      curr <- plot_df[i, ]
      next_row <- plot_df[i + 1, ]
      
      mid <- (curr$end + next_row$start) / 2
      
      # Add horizontal line (from start to midpoint)
      step_df <- rbind(step_df, data.frame(x = curr$start, y = curr$read_ratio))
      step_df <- rbind(step_df, data.frame(x = mid, y = curr$read_ratio))
      
      # Add vertical transition
      step_df <- rbind(step_df, data.frame(x = mid, y = next_row$read_ratio))
    }
    
    # Add the final horizontal segment to the last region's end
    step_df <- rbind(step_df,
                     data.frame(x = plot_df$start[nrow(plot_df)], y = plot_df$read_ratio[nrow(plot_df)]),
                     data.frame(x = plot_df$end[nrow(plot_df)], y = plot_df$read_ratio[nrow(plot_df)]))
    
    read_ratio_plot <- ggplot(plot_df) +
      # Original horizontal bars
      geom_segment(aes(x = start, xend = end, y = read_ratio, yend = read_ratio),
                   color = "#1f77b4", size = 1.2) +
      
      # 90-degree angle line path
      geom_path(data = step_df, aes(x = x, y = y), color = "#1f77b4", size = 0.4) +
      
      # Reference line
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      
      # Duplication and deletion threshold lines
      geom_hline(yintercept = 1.5, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
      
      # Labels for the duplication and deletion lines
      annotate("text", x = user_start, y = 1.5, label = "Read Ratio = 1.5 (Duplication)",
               hjust = 0, vjust = -0.5, color = "red", size = 3.5) +
      annotate("text", x = user_start, y = 0.5, label = "Read Ratio = 0.5 (Deletion)",
               hjust = 0, vjust = -0.5, color = "red", size = 3.5) +
      # Translucent box spanning user_start to user_end
      geom_rect(aes(xmin = cnv_start, xmax = cnv_end, ymin = -Inf, ymax = Inf),
                fill = "blue", alpha = 0.01) +
      
      theme_bw() +
      labs(
        title = paste0("NGS Read Ratio and aCGH Log2Ratio for Sample ", str_extract(target_sample_dots, "DNS\\d+")),
        subtitle = paste0("chr", user_chr, ":", user_start, "-", user_end),
        x = NULL,
        y = "Read Ratio"
      ) +
      scale_x_continuous(labels = scales::comma) +
      scale_y_continuous(
        breaks = c(2.0, 1.75, 1.5, 1.25, 1.0, 0.75, 0.5, 0.25, 0),
        limits = c(0, 2.3)
      ) +
      xlim(user_start, user_end) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
      )
  } 
  
  if ("exon_plot" %in% plot_types) { 
    print(paste0("plotting exon track for ", target_sample, "_", gene_symbol))
    
    # --- Parameters ---
    clinically_relevant_transcripts = list(BRCA1 = "NM_007294.4",
                                           BRCA2 = "NM_000059.4",
                                           PALB2 = "NM_024675.4",
                                           MSH2 = "NM_000251.3",
                                           MSH6 = "NM_000179.3",
                                           MLH1 = "NM_000249.4",
                                           APOB = "NM_000384.3",
                                           PCSK9 = "NM_174936.4",
                                           LDLR = "NM_000527.5")
    
    
    refseq_id <- clinically_relevant_transcripts[[gene_symbol]]
    
    # Load cached exon data
    exons_combined <- readRDS("/path/to/project/Aim_2/CNV-plot/cached_exon_data.rds")
    
    # Filter for gene of interest and region
    exons_region <- exons_combined %>%
      filter(gene == gene_symbol) %>%
      filter(chromosome_name == as.character(user_chr)) %>%
      filter(exon_chrom_start >= user_start & exon_chrom_end <= user_end)
    
    # 1. Sort exons by start position
    exons_sorted <- exons_region %>% arrange(exon_chrom_start)
    
    # 2. Determine strand and arrow symbol
    strand <- unique(exons_sorted$strand)
    arrow_symbol <- ifelse(strand == -1, "<", ">")
    
    # 3. Set spacing and y position
    arrow_spacing <- 500  # You can tweak this
    arrow_y <- 0.5
    
    # 4. Build arrow data frame between exons
    arrow_df <- data.frame()
    
    for (i in 1:(nrow(exons_sorted) - 1)) {
      intron_start <- exons_sorted$exon_chrom_end[i] + 1
      intron_end   <- exons_sorted$exon_chrom_start[i + 1] - 1
      
      if (intron_start < intron_end) {
        arrow_positions <- seq(intron_start, intron_end, by = arrow_spacing)
        arrow_df <- rbind(
          arrow_df,
          data.frame(x = arrow_positions, y = arrow_y, label = arrow_symbol)
        )
      }
    }
    
    # Define transcript_id for the current gene (e.g., "NM_000384.3")
    transcript_id <- clinically_relevant_transcripts[[gene_symbol]]
    
    # Compute midpoint of the gene region for positioning the label
    gene_mid <- exons_region %>%
      summarize(mid = mean(range(exon_chrom_start, exon_chrom_end))) %>%
      pull(mid)
    
    # Plot with arrows and label
    exon_plot <- ggplot(exons_region) +
      geom_hline(yintercept = arrow_y, color = "black", size = 0.4) +
      
      # Arrow path between exons
      geom_text(data = arrow_df,
                aes(x = x, y = y, label = label),
                size = 3, color = "black", fontface = "bold") +
      
      # Exon boxes on top
      geom_rect(aes(xmin = exon_chrom_start, xmax = exon_chrom_end,
                    ymin = 0.4, ymax = 0.6, fill = gene),
                color = "black", alpha = 1) +
      
      # Transcript label
      geom_text(data = data.frame(x = gene_mid, y = 0.7,
                                  label = paste0(gene_symbol, " (", transcript_id, ")")),
                aes(x = x, y = y, label = label),
                size = 3, fontface = "bold") +
      
      theme_void() +
      scale_fill_brewer(palette = "Set2") +
      theme(
        legend.position = "none",
        plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
      ) +
      xlim(user_start, user_end)
  }
  if ("raw_depth_plot" %in% plot_types) { 
    print(paste0("plotting raw depth track for ", target_sample, "_", gene_symbol))
    
    roi <- GRanges(
      seqnames = as.character(paste0("chr", user_chr)),
      ranges = IRanges(start = user_start, end = user_end)
    )
    
    # Full path to BAM and index
    val_bam_dir <- "/path/to/project/Aim_1/Validations/bam_and_vcf"
    bam_path <- file.path(val_bam_dir, gsub("\\.bam$", ".bam", target_sample))
    bam_index <- paste0(bam_path, ".bai")
    
    # Create BamFile object
    bam_file <- BamFile(bam_path, index = bam_index)
    
    # Extract reads overlapping the region of interest
    reads <- readGAlignments(bam_file, param = ScanBamParam(which = roi))
    
    # Compute coverage
    cov <- coverage(reads)
    
    # Subset to your region of interest
    cov_region <- cov[[paste0("chr", user_chr)]][user_start:user_end]
    
    # Convert to data frame
    raw_depth_df <- data.frame(
      position = seq(user_start, user_end),
      depth = as.numeric(cov_region)
    )
    
    raw_depth_plot <- ggplot(raw_depth_df, aes(x = position, y = depth)) +
      geom_line(color = "#009E73") +
      labs(
        title = "Raw Read Depth",
        y = "Coverage",
        x = NULL
      ) +
      theme_bw() +
      scale_x_continuous(labels = scales::comma) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
      )
  }
  
  ### Dan's request ####
  
  coords <- gsub("^(\\d+):(\\d+)-(\\d+)$", "chr\\1_\\2_\\3", cnv_row$coordinates)
  
  # Combine vertically
  combined_plot_withRD <- raw_depth_plot / read_ratio_plot  / exon_plot +
    plot_layout(heights = c(2, 5.8, 0.8))
  
  # Use save_folder in ggsave filename
  ggsave(
    filename = paste0(save_folder, str_extract(target_sample, "DNS\\d+"), "_", gene_symbol, "_", coords, "_read_ratio_with_rawRD.jpeg"),
    plot = combined_plot_withRD,
    width = 12,
    height = 6,
    units = "in"
  )
  
}



# # Loop through rows in vscnv_reorder_perrun and call the function
# for (i in seq_len(nrow(validated_res))) {
#   row <- validated_res[i, ]  # This keeps the original data types
#   
#   # Convert the one-row data.frame to a named list for your function (if needed)
#   row_list <- as.list(row)
#   
#   plot_cnv_for_sample(row_list, read_count_dir, bed_name, bed_file, 
#                       plot_types = c("read_ratio_plot", "exon_plot", "raw_depth_plot"))
# }


for (i in seq_len(nrow(validated_res_clean))) {
  row <- validated_res_clean[i, ]
  print(row$Sample)
  row_list <- as.list(row)
  
  # Determine folder based on Verified_20250730 column
  if (!is.null(row$Verified_20250730) && row$Verified_20250730 == "YES") {
    folder <- "/path/to/project/Aim_2/CNV-plot/concordant/"
  } else {
    folder <- "/path/to/project/Aim_2/CNV-plot/discordant/"
  }
  
  # Create folder if doesn't exist
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  
  plot_cnv_for_sample(row_list, read_count_dir, bed_name, bed_file, 
                      plot_types = c("read_ratio_plot", "exon_plot", "raw_depth_plot"), gc_correction = FALSE,
                      save_folder = folder)
}


