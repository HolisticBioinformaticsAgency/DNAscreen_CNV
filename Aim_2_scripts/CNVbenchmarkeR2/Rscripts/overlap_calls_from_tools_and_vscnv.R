library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tibble)
library(scales)
library(plotly)
library(rlang)
library(GenomicRanges)

dir <- "/path/to/project/Aim_1/CNV_calling/cnvbenchmarker/CNVbenchmarkeR2_test/output"
directory_bam <- "/path/to/project/Aim_1/Validations/bam_and_vcf"
directory_vcf <- "/path/to/project/Aim_1/Validations/bam_and_vcf"

# dir <- "/path/to/project/Aim_1/CNV_calling/cnvbenchmarker/CNVbenchmarkeR2_test/output_run10"
# directory_bam <- "/path/to/project/Aim_1/CNV_calling/exome_depth/bams_run31"
# directory_vcf <- "/path/to/project/Aim_1/CNV_calling/VarSeq_CNV/vcfs_run31"

# List all files in the directory with the extension '*.bam'
bam_files <- list.files(directory_bam, pattern = ".bam$", full.names = TRUE)
bai_files <- list.files(directory_bam, pattern = ".bai$", full.names = TRUE)
# vcf files are imported since only samples with both bam and vcf passed QC
vcf_files <- list.files(directory_vcf, pattern = ".vcf.gz$", full.names = TRUE)
vcf_filenames <- basename(vcf_files)

extracted <- sub("\\.sorted\\.vcf\\.gz$", "", vcf_filenames)

# Extract sample IDs from bam_files using regex
sample_ids <- sub(".*/([^/]+)\\.hq\\.sorted\\.marked\\.bam", "\\1", bam_files)
sample_bai_ids <- sub(".*/([^/]+)\\.hq\\.sorted\\.marked\\.bam.bai", "\\1", bai_files)

# Check if sample IDs exist in 'extracted' and 'sample_bai_ids'
passqc_samples <- sample_ids %in% extracted & sample_ids %in% sample_bai_ids

# Filter sample_ids based on samples that passed QC
sample_ids <- sample_ids[passqc_samples]


atlascnv <- read.csv(paste0(dir, "/atlasCNV-dataset1/cnv_calls.tsv"), sep = "\t")
clearcnv <- read.csv(paste0(dir, "/clearCNV-dataset1/cnv_calls.tsv"), sep = "\t")
clincnv <- read.csv(paste0(dir, "/clincnv-dataset1/cnvs_summary.tsv"), sep = "\t")
cnvkit <- read.csv(paste0(dir, "/cnvkit-dataset1/cnvs_summary.tsv"), sep = "\t")
cobalt <- read.csv(paste0(dir, "/cobalt-dataset1/cnvs_summary.tsv"), sep = "\t")
codex2 <- read.csv(paste0(dir, "/codex2-dataset1/cnvFounds.csv"), sep = "\t")
decon <- read.csv(paste0(dir, "/decon-dataset1/calls_all.txt"), sep = "\t")
exomedepth <- read.csv(paste0(dir, "/exomedepth-dataset1/all_cnv_calls.txt"), sep = "\t")
gatk <- read.csv(paste0(dir, "/germlineCNVcaller-dataset1/all_cnv_calls.txt"), sep = "\t")
panelcn <- read.csv(paste0(dir, "/panelcnmops-dataset1/cnvFounds.txt"), sep = "\t") #does not have a quality metric
viscap <- read.csv(paste0(dir, "/viscap-dataset1/cnvs_summary.tsv"), sep = "\t")


### This part is mainly for GATK since it doesn't assign gene names
### and for acgh (Dan gave coordinates different from the varseq calls so i can't use the gene column from varseq)
# Read the gene coordinates
gene_data <- read.table("/path/to/project/Aim_1/CNV_calling/CI_plots/dnascreen_genes.txt", header = FALSE, col.names = c("chr", "start", "end", "Gene.Names"))
gene_data$chr <- str_replace(gene_data$chr, "^chr", "")
# Convert gene_data to GRanges
gene_ranges <- GRanges(
  seqnames = gene_data$chr,
  ranges = IRanges(start = gene_data$start, end = gene_data$end),
  gene_name = gene_data$Gene.Names
)

### Include bed coordinates of round one validated CNVs
# Read validated CNV data (that are CONFIRMED (no maybe or likely)) and clean
validated_res <- read.csv("/path/to/project/Aim_1/Validations/acgh_validation_status_files/acgh_round_one_confirmed_cnvs_dan.csv")

validated_res$Man..cur..Chr <- gsub("^chr", "", validated_res$Man..cur..Chr)
validated_res$Man..cur..Chr[validated_res$Man..cur..Chr == ""] <- NA
validated_res[] <- lapply(validated_res, function(x) {
  if (is.character(x)) x[x == ""] <- NA
  return(x)
})

# Filter and clean subset of data
acgh <- validated_res %>%
  filter(
    Verified_20250730 == "YES"
      # (!is.na(Manual.curation.Start) & !is.na(Manual.curation.Stop) & !is.na(Man..cur..Chr)) # Dan said concentrate on YES only
  ) %>%
  mutate(
    verified_by_dan = Verified_20250730,
    cn = recode(cn, "1" = "deletion", "3" = "duplication")
  ) %>%
  select(sample, Gene.Names, verified_by_dan, Manual.curation.Start, Manual.curation.Stop, Man..cur..Chr, cn)

# Convert manual curation coordinates to GRanges
acgh_gr <- GRanges(
  seqnames = acgh$Man..cur..Chr,
  ranges = IRanges(start = acgh$Manual.curation.Start, end = acgh$Manual.curation.Stop)
)

# Perform overlaps
overlaps <- findOverlaps(acgh_gr, gene_ranges)

# Extract overlapping gene names
acgh$Gene.from.coords <- NA  # placeholder
acgh$Gene.from.coords[queryHits(overlaps)] <- gene_ranges$gene_name[subjectHits(overlaps)]


# Load per-run 1392 vs-cnv calls
vs_cnv_calls_perrun <- read.csv("/path/to/project/Aim_1/CNV_calling/CI_plots/combined_vs_tables/varseq_sample_cnv_calls_combined_all_runs_v3_with_dan_verification_cleaned.csv")
# Load validation set only vs-cnv calls
vs_cnv_calls_val <- read.csv("/path/to/project/Aim_1/CNV_calling/CI_plots/combined_vs_tables/varseq_sample_cnv_calls_combined_with_validation_samples_as_ref_set_v2.csv")

# harmonise all tool output tables
safe_na_num <- function() as.numeric(NA)
safe_na_character <- function() as.character(NA)

# VS-CNV perrun
vscnv_reorder_perrun <- vs_cnv_calls_perrun %>%
  dplyr::rename(Sample = sample,
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
  # also need to limit to just samples that are in the validation bam folder!
  filter(Sample %in% sample_ids) %>%
  filter(!is.na(Gene)) %>% ## Filter out any CNVs in ATM ROI
  select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

  
  
vscnv_harmonised_perrun <- vscnv_reorder_perrun %>%
  mutate(Tool = "VS-CNV_per_run", read_ratio = round(Avg.Read.Ratio.of.CNV, 2), 
         Chromosome = as.character(Chromosome),
         log2ratio = safe_na_num(), 
         quality_score = paste("(p-value)", adjusted_p_value)) %>%
  select(Tool, Sample, Chromosome, Start, End, Gene, 
         CNV_Type, read_ratio, log2ratio, quality_score)

# VS-CNV validation only as ref set
vscnv_reorder_val <- vs_cnv_calls_val %>%
  dplyr::rename(
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
  # also need to limit to just samples that are in the bam folder
  filter(Sample %in% sample_ids) %>%
  filter(!is.na(Gene)) %>%
  select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())


vscnv_harmonised_val <- vscnv_reorder_val %>%
  mutate(Tool = "VS-CNV_val", read_ratio = round(Avg.Read.Ratio.of.CNV, 2), 
         Chromosome = as.character(Chromosome),
         log2ratio = safe_na_num(), 
         quality_score = paste("(p-value)", adjusted_p_value)) %>%
  select(Tool, Sample, Chromosome, Start, End, Gene, 
         CNV_Type, read_ratio, log2ratio, quality_score)


atlascnv_reorder <- atlascnv %>%
  dplyr::rename(
    Sample = sample,
    Chromosome = chr,
    Start = start,
    End = end,
    Gene = RefSeq,
    CNV_Type = CNV.type
  ) %>%
  arrange(desc(abs(c_Score))) %>% 
  mutate(Sample = sub("\\.hq\\.sorted\\.marked.*", "", Sample),
         Chromosome = gsub("^chr", "", Chromosome)) %>%
  # also need to limit to just samples that are in the bam folder
  filter(Sample %in% sample_ids) %>%
  filter(!is.na(Gene)) %>%
  select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

write.table(atlascnv_reorder, paste0(dir, "/atlasCNV-dataset1/atlascnv_reorder.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

atlascnv_harmonised <- atlascnv_reorder %>%
  mutate(Tool = "atlasCNV", read_ratio = safe_na_num(), 
         log2ratio = round(log2R, 2), quality_score = paste("(c_Score)", round(c_Score, 2))) %>%
  select(Tool, Sample, Chromosome, Start, End, Gene, 
         CNV_Type, read_ratio, log2ratio, quality_score)

clearcnv_reorder <- clearcnv %>%
  dplyr::rename(
    Sample = sample,
    Chromosome = chr,
    Start = start,
    End = end,
    Gene = gene,
    CNV_Type = CNV.type
  ) %>%
  arrange(desc(score)) %>% 
  mutate(Sample = sub("\\.hq\\.sorted\\.marked.*", "", Sample),
         Chromosome = gsub("^chr", "", Chromosome)) %>%
  # also need to limit to just samples that are in the bam folder
  filter(Sample %in% sample_ids) %>%
  filter(!is.na(Gene)) %>%
  select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

write.table(clearcnv_reorder, paste0(dir, "/clearCNV-dataset1/clearcnv_reorder.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

clearcnv_harmonised <- clearcnv_reorder %>%
  mutate(Tool = "clearCNV", read_ratio = round(ratio, 2), 
         log2ratio = safe_na_num(), quality_score = paste("(score)", round(score, 2))) %>%
  select(Tool, Sample, Chromosome, Start, End, Gene, 
         CNV_Type, read_ratio, log2ratio, quality_score)

clincnv_reorder <- clincnv %>%
  dplyr::rename(
    Sample = Sample,
    Chromosome = chr,
    Start = start,
    End = end,
    Gene = genes,
    CNV_Type = CNV.type
  ) %>%
  arrange(desc(loglikelihood)) %>%
  mutate(Sample = sub("\\.hq\\.sorted\\.marked.*", "", Sample),
         Chromosome = gsub("^chr", "", Chromosome)) %>%
  # also need to limit to just samples that are in the bam folder
  filter(Sample %in% sample_ids) %>%
  filter(!is.na(Gene)) %>%
  select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

write.table(clincnv_reorder, paste0(dir, "/clincnv-dataset1/clincnv_reorder.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

clincnv_harmonised <- clincnv_reorder %>%
  mutate(Tool = "clinCNV", read_ratio = safe_na_num(), 
         log2ratio = safe_na_num(), 
         quality_score = paste("(loglikelihood)", loglikelihood)) %>%
  select(Tool, Sample, Chromosome, Start, End, Gene, 
         CNV_Type, read_ratio, log2ratio, quality_score)

cnvkit_reorder <- cnvkit %>%
  dplyr::rename(
    Sample = Sample,
    Chromosome = seqnames,
    Start = start_gene,
    End = end_gene,
    Gene = name,
    CNV_Type = CNV.type
  ) %>%
  arrange(desc(log2)) %>%
  mutate(Sample = sub("\\.hq\\.sorted\\.marked.*", "", Sample),
         Chromosome = gsub("^chr", "", Chromosome)) %>%
  # also need to limit to just samples that are in the bam folder
  filter(Sample %in% sample_ids) %>%
  filter(!is.na(Gene)) %>%
  select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

write.table(cnvkit_reorder, paste0(dir, "/cnvkit-dataset1/cnvkit_reorder.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

cnvkit_harmonised <- cnvkit_reorder %>%
  mutate(Tool = "cnvkit", read_ratio = safe_na_num(), 
         log2ratio = round(log2 ,2), quality_score = safe_na_character()) %>%
  select(Tool, Sample, Chromosome, Start, End, Gene, 
         CNV_Type, read_ratio, log2ratio, quality_score)

cobalt_reorder <- cobalt %>%
  dplyr::rename(
    Sample = Sample,
    Chromosome = seqnames,
    Start = start_gene,
    End = end_gene,
    Gene = name,
    CNV_Type = CNV.type
  ) %>%
  arrange(desc(quality)) %>%
  mutate(Sample = sub("\\.hq\\.sorted\\.marked.*", "", Sample),
         Chromosome = gsub("^chr", "", Chromosome)) %>%
  # also need to limit to just samples that are in the bam folder
  filter(Sample %in% sample_ids) %>%
  filter(!is.na(Gene)) %>%
  select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

write.table(cobalt_reorder, paste0(dir, "/cobalt-dataset1/cobalt_reorder.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

cobalt_harmonised <- cobalt_reorder %>%
  mutate(Tool = "cobalt", read_ratio = safe_na_num(), 
         log2ratio = safe_na_num(), 
         quality_score = paste("(quality)", quality)) %>%
  select(Tool, Sample, Chromosome, Start, End, Gene, 
         CNV_Type, read_ratio, log2ratio, quality_score)

codex2_reorder <- codex2 %>%
  dplyr::rename(
    Sample = sample_name,
    Chromosome = chr,
    Start = st_bp,
    End = ed_bp,
    Gene = gene,
    CNV_Type = cnv
  ) %>%
  arrange(desc(mBIC)) %>%
  mutate(Sample = sub("\\.hq\\.sorted\\.marked.*", "", Sample),
         Chromosome = gsub("^chr", "", Chromosome)) %>%
  # also need to limit to just samples that are in the bam folder
  filter(Sample %in% sample_ids) %>%
  filter(!is.na(Gene)) %>%
  select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

write.table(codex2_reorder, paste0(dir, "/codex2-dataset1/codex2_reorder.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

codex2_harmonised <- codex2_reorder %>%
  mutate(Tool = "codex2", read_ratio = round(copy_no/2, 2), 
         log2ratio = safe_na_num(), quality_score = paste("(mBIC)", mBIC)) %>%
  select(Tool, Sample, Chromosome, Start, End, Gene, 
         CNV_Type, read_ratio, log2ratio, quality_score)

decon_reorder <- decon %>%
  dplyr::rename(
    Sample = Sample,
    Chromosome = Chromosome,
    Start = Start,
    End = End,
    Gene = Gene,
    CNV_Type = CNV.type
  ) %>%
  arrange(desc(BF)) %>%
  mutate(Sample = sub("\\.hq\\.sorted\\.marked.*", "", Sample),
         Chromosome = gsub("^chr", "", Chromosome)) %>%
  # also need to limit to just samples that are in the bam folder
  filter(Sample %in% sample_ids) %>%
  filter(!is.na(Gene)) %>%
  select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

write.table(decon_reorder, paste0(dir, "/decon-dataset1/decon_reorder.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

decon_harmonised <- decon_reorder %>%
  mutate(Tool = "decon", read_ratio = Reads.ratio, 
         log2ratio = safe_na_num(), quality_score = paste("(BF)", BF)) %>%
  select(Tool, Sample, Chromosome, Start, End, Gene, 
         CNV_Type, read_ratio, log2ratio, quality_score)

exomedepth_reorder <- exomedepth %>%
  dplyr::rename(
    Sample = sample,
    Chromosome = chromosome,
    Start = start,
    End = end,
    Gene = Gene,
    CNV_Type = type
  ) %>%
  arrange(desc(BF)) %>%
  mutate(Sample = sub("\\.hq\\.sorted\\.marked.*", "", Sample),
         Chromosome = gsub("^chr", "", Chromosome)) %>%
  # also need to limit to just samples that are in the bam folder
  filter(Sample %in% sample_ids) %>%
  filter(!is.na(Gene)) %>%
  select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

write.table(exomedepth_reorder, paste0(dir, "/exomedepth-dataset1/exomedepth_reorder.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

exomedepth_harmonised <- exomedepth_reorder %>%
  mutate(Tool = "exomedepth", read_ratio = reads.ratio, 
         log2ratio = safe_na_num(), quality_score = paste("(BF)", BF)) %>%
  select(Tool, Sample, Chromosome, Start, End, Gene, 
         CNV_Type, read_ratio, log2ratio, quality_score)

gatk_reorder <- gatk %>%
  mutate(Gene = safe_na_character()) %>%
  dplyr::rename(
    Sample = sample,
    Chromosome = CHROM,
    Start = POS,
    End = END,
    Gene = Gene,
    CNV_Type = copy_number
  ) %>%
  arrange(desc(gt_QS)) %>%
  mutate(Sample = sub("\\.hq\\.sorted\\.marked.*", "", Sample),
         Chromosome = gsub("^chr", "", Chromosome)) %>%
  # also need to limit to just samples that are in the bam folder
  filter(Sample %in% sample_ids) %>%
  select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

# Convert gatk_reorder to GRanges
gatk_ranges <- GRanges(
  seqnames = gatk_reorder$Chromosome,
  ranges = IRanges(start = gatk_reorder$Start, end = gatk_reorder$End)
)

# Find overlaps
hits <- findOverlaps(gatk_ranges, gene_ranges)

# Assign gene names to gatk_reorder
gatk_reorder$Gene <- NA_character_
gatk_reorder$Gene[queryHits(hits)] <- gene_ranges$gene_name[subjectHits(hits)]
gatk_reorder <- gatk_reorder %>%
  filter(!is.na(Gene))

write.table(gatk_reorder, paste0(dir, "/germlineCNVcaller-dataset1/gatk_reorder.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

gatk_harmonised <- gatk_reorder %>%
  mutate(Tool = "gatk", read_ratio = safe_na_num(), 
         log2ratio = safe_na_num(), quality_score = paste("(gt_QS)", gt_QS)) %>%
  select(Tool, Sample, Chromosome, Start, End, Gene, 
         CNV_Type, read_ratio, log2ratio, quality_score)

panelcn_reorder <- panelcn %>%
  dplyr::rename(
    Sample = Sample,
    Chromosome = Chr,
    Start = Start,
    End = End,
    Gene = Gene,
    CNV_Type = CNV.type
  ) %>%
  mutate(Sample = sub("\\.hq\\.sorted\\.marked.*", "", Sample),
         Chromosome = gsub("^chr", "", Chromosome)) %>%
  # also need to limit to just samples that are in the bam folder
  filter(Sample %in% sample_ids) %>%
  filter(!is.na(Gene)) %>%
  select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

write.table(panelcn_reorder, paste0(dir, "/panelcnmops-dataset1/panelcn_reorder.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

panelcn_harmonised <- panelcn_reorder %>%
  mutate(Tool = "panelcn.MOPS", read_ratio = safe_na_num(), 
         log2ratio = safe_na_num(), quality_score = safe_na_character()) %>%
  select(Tool, Sample, Chromosome, Start, End, Gene, 
         CNV_Type, read_ratio, log2ratio, quality_score)

viscap_reorder <- viscap %>%
  dplyr::rename(
    Sample = Sample,
    Chromosome = seqnames,
    Start = start_gene,
    End = end_gene,
    Gene = name,
    CNV_Type = copy_number
  ) %>%
  arrange(desc(abs(quality))) %>% # quality here is actually log2ratio
  mutate(Sample = sub("\\.hq\\.sorted\\.marked.*", "", Sample),
         Chromosome = gsub("^chr", "", Chromosome)) %>%
  # also need to limit to just samples that are in the bam folder
  filter(Sample %in% sample_ids) %>%
  filter(!is.na(Gene)) %>%
  select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

write.table(viscap_reorder, paste0(dir, "/viscap-dataset1/viscap_reorder.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

viscap_harmonised <- viscap_reorder %>%
  mutate(Tool = "VisCap", read_ratio = safe_na_num(), 
         log2ratio = quality, quality_score = safe_na_character()) %>%
  select(Tool, Sample, Chromosome, Start, End, Gene, 
         CNV_Type, read_ratio, log2ratio, quality_score)

## ROUND ONE ACGH VERIFIED == YES BY DAN! ###

acgh_reorder <- acgh %>%
  dplyr::rename(
    Sample = sample,
    Chromosome = Man..cur..Chr,
    Start = Manual.curation.Start,
    End = Manual.curation.Stop,
    Gene = Gene.from.coords,
    CNV_Type = cn
  ) %>%
  mutate(Sample = sub("\\.hq\\.sorted\\.marked.*", "", Sample),
         Chromosome = gsub("^chr", "", Chromosome)) %>%
  # also need to limit to just samples that are in the bam folder
  filter(Sample %in% sample_ids) %>%
  filter(!is.na(Gene)) %>%
  select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

write.table(acgh_reorder, paste0(dir, "/acgh_round_one/acgh_reorder.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

acgh_harmonised <- acgh_reorder %>%
  mutate(Tool = "acgh", read_ratio = safe_na_num(), 
         log2ratio = safe_na_num(), quality_score = safe_na_character()) %>%
  select(Tool, Sample, Chromosome, Start, End, Gene, 
         CNV_Type, read_ratio, log2ratio, quality_score)



############ Code below writes a bed file of all the cnv calls on a per-sample basis ##############

for (sample_id in sample_ids) {
  # Filter harmonised tables for this sample
  vs_cnv_row_perrun <- vscnv_harmonised_perrun %>% filter(Sample == sample_id)
  vs_cnv_row_val    <- vscnv_harmonised_val %>% filter(Sample == sample_id)
  atlascnv_row      <- atlascnv_harmonised %>% filter(Sample == sample_id)
  clearcnv_row      <- clearcnv_harmonised %>% filter(Sample == sample_id)
  clincnv_row       <- clincnv_harmonised %>% filter(Sample == sample_id)
  cnvkit_row        <- cnvkit_harmonised %>% filter(Sample == sample_id)
  cobalt_row        <- cobalt_harmonised %>% filter(Sample == sample_id)
  codex2_row        <- codex2_harmonised %>% filter(Sample == sample_id)
  decon_row         <- decon_harmonised %>% filter(Sample == sample_id)
  exomedepth_row    <- exomedepth_harmonised %>% filter(Sample == sample_id)
  gatk_row          <- gatk_harmonised %>% filter(Sample == sample_id)
  panelcn_row       <- panelcn_harmonised %>% filter(Sample == sample_id)
  viscap_row        <- viscap_harmonised %>% filter(Sample == sample_id)
  acgh_row          <- acgh_harmonised %>% filter(Sample == sample_id)
  
  # Combine into one data frame
  overlap <- bind_rows(
    vs_cnv_row_perrun, vs_cnv_row_val, atlascnv_row, clearcnv_row, clincnv_row, cnvkit_row, cobalt_row,
    codex2_row, decon_row, exomedepth_row, gatk_row, panelcn_row, viscap_row, acgh_row
  )
  
  # Define output directory for this sample
  output_sample_dir <- file.path(dir, "per_sample_cnv_calls", sample_id)
  if (!dir.exists(output_sample_dir)) {
    dir.create(output_sample_dir, recursive = TRUE)
  }
  
  # Write combined table
  write.table(
    overlap,
    file = file.path(output_sample_dir, "cnv_calls.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  # Optionally, write BED files too
  write_bed_per_tool(overlap, output_dir = file.path(output_sample_dir, "bed_files"))
}


# List of harmonised data frames and their quality score variables
### to be used in tool_roi_gene_accuracy.R, jitter_plot_hgsa and which_tool_calls_cnv.R
harmonised_list <- list(
  "VS-CNV_perrun" = list(data = vscnv_reorder_perrun, quality_col = "log_adjusted_p_value"),
  "VS-CNV_val" = list(data = vscnv_reorder_val, quality_col = "log_adjusted_p_value"),
  "atlasCNV"   = list(data = atlascnv_reorder,     quality_col = "c_Score"),
  "clearCNV"   = list(data = clearcnv_reorder,     quality_col = "score"),
  "clinCNV"    = list(data = clincnv_reorder,      quality_col = "loglikelihood"),
  "cnvkit"     = list(data = cnvkit_reorder,       quality_col = "log2"),
  "cobalt"     = list(data = cobalt_reorder,       quality_col = "quality"),
  "codex2"     = list(data = codex2_reorder,       quality_col = "mBIC"),
  "decon"      = list(data = decon_reorder,        quality_col = "BF"),
  "exomedepth" = list(data = exomedepth_reorder,   quality_col = "BF"),
  "gatk"       = list(data = gatk_reorder,         quality_col = "gt_QS"),
  "VisCap"     = list(data = viscap_reorder,       quality_col = "quality"),
  "panelcn"     = list(data = panelcn_reorder,     quality_col = "None")
)
# Note: panelcn does not output any quality scores...


