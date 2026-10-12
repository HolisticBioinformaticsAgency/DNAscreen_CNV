library(dplyr)
library(readr)
library(stringr)

# Function to assign gene names to CNV calls
assign_gene_names <- function(cnv_calls, gene_data) {
  # Create GenomicRanges for CNV calls
  cnv_gr <- GRanges(
    seqnames = cnv_calls$VS_Call_chr,
    ranges = IRanges(
      start = cnv_calls$VS_Call_start,
      end = cnv_calls$VS_Call_end
    )
  )
  
  # Create GenomicRanges for gene data
  gene_gr <- GRanges(
    seqnames = gene_data$chr,
    ranges = IRanges(start = gene_data$start, end = gene_data$end),
    Gene.Names = gene_data$Gene.Names
  )
  
  # Find overlaps between CNV regions and gene regions
  overlaps <- findOverlaps(cnv_gr, gene_gr)
  
  # Initialize the 'Gene.Names' column with NA
  cnv_calls$Gene.Names <- NA
  
  # Assign gene names to the overlapping CNV regions
  cnv_calls$Gene.Names[queryHits(overlaps)] <- gene_gr$Gene.Names[subjectHits(overlaps)]
  
  return(cnv_calls)
}

write_bed_per_tool <- function(cnv_data, output_dir = "bed_files") {
  library(dplyr)
  library(stringr)
  library(readr)
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Add 'chr' prefix to Chromosome
  cnv_data <- cnv_data %>%
    mutate(chr = paste0("chr", Chromosome))
  
  # Label the ratio or log2ratio explicitly
  cnv_data <- cnv_data %>%
    mutate(
      ratio_label = case_when(
        !is.na(read_ratio) ~ paste0("ratio:", read_ratio),
        !is.na(log2ratio) ~ paste0("log2ratio:", log2ratio),
        TRUE ~ NA_character_
      )
    )
  
  # Parse quality score into human-readable labels
  cnv_data <- cnv_data %>%
    mutate(
      quality_label = case_when(
        str_detect(quality_score, "p-value") ~ paste0("pval:", str_extract(quality_score, "(?<=\\))\\s*.*$")),
        str_detect(quality_score, "c_Score") ~ paste0("cScore:", str_extract(quality_score, "-?[0-9.]+")),
        str_detect(quality_score, "score") & !str_detect(quality_score, "c_Score") ~ paste0("score:", str_extract(quality_score, "[0-9.]+")),
        str_detect(quality_score, "loglikelihood") ~ paste0("LL:", str_extract(quality_score, "[0-9.]+")),
        str_detect(quality_score, "quality") ~ paste0("qual:", str_extract(quality_score, "[0-9.]+")),
        str_detect(quality_score, "mBIC") ~ paste0("mBIC:", str_extract(quality_score, "[0-9.]+")),
        str_detect(quality_score, "BF") ~ paste0("BF:", str_extract(quality_score, "[0-9.]+")),
        str_detect(quality_score, "gt_QS") ~ paste0("gt_QS:", str_extract(quality_score, "[0-9.]+")),
        TRUE ~ NA_character_
      )
    )
  
  # Combine into 'name' column for BED format
  cnv_data <- cnv_data %>%
    mutate(
      name = str_c(
        Start, "-", End,
        ifelse(!is.na(CNV_Type), paste0("_", CNV_Type), ""),
        ifelse(!is.na(ratio_label), paste0("_", ratio_label), ""),
        ifelse(!is.na(quality_label), paste0("_", quality_label), "")
      )
    )
  
  # Write BED files per Tool
  tools <- unique(cnv_data$Tool)
  
  for (tool in tools) {
    tool_df <- cnv_data %>%
      filter(Tool == tool) %>%
      select(chr, Start, End, name) %>%
      arrange(Start)
    
    # Output path
    output_path <- file.path(output_dir, paste0(tool, ".bed"))
    
    # Write BED file
    write_tsv(tool_df, file = output_path, col_names = FALSE)
  }
}
