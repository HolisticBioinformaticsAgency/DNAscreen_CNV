library(dplyr)
library(tidyr)

# # --- 1. Read in validated TP rows ---
# validated_run <- read.csv(
#   '/path/to/project/Aim_1/CNV_calling/cnvbenchmarker/CNVbenchmarkeR2_test/examples/validated_results_val_samples_ICR_format.tsv',
#   sep = '\t'
# )
# tp_rows <- validated_run  # assuming it contains all validated TPs

# --- 2. Get gene info ---
dnascreen_genes <- read.csv(
  "/path/to/project/Aim_1/CNV_calling/cnvbenchmarker/CNVbenchmarkeR2_test/bed_file/dnascreen_genes_edited_to_match_covered.bed",
  sep = "\t", header = FALSE
)
colnames(dnascreen_genes) <- c("chr", "start", "end", "gene", "length")
dnascreen_genes$chr <- as.integer(gsub("chr", "", dnascreen_genes$chr))

# --- 3. Get sample IDs from samples that have already been validated ---
# bam_dir <- "/path/to/project/Aim_1/Validations/bam_and_vcf"
# # bam_dir <- "/path/to/project/Aim_1/CNV_calling/exome_depth/bams_run31"
# bam_files <- list.files(path = bam_dir, pattern = "\\.bam$", full.names = FALSE)
# sample_ids <- sub("\\.bam$", "", bam_files)

### Now using validated_res
validated_res <- read.csv("/path/to/project/Aim_1/Validations/acgh_validation_status_files/acgh_round_one_confirmed_cnvs_dan.csv")
sample_ids <- unique(validated_res$sample)
validated_res <- validated_res %>%
  mutate(Man..cur..Chr = as.integer(gsub("^chr", "", Man..cur..Chr)))

validated_res <- validated_res %>%
  mutate(Type.of.CNV = recode(Type.of.CNV,
                              "Het Deletion" = "Deletion",
                              "Duplicate" = "Duplication"))

tp_rows <- validated_res %>%
  filter(Verified_20250730 == "YES") %>%
  mutate(SampleID = sample,
         ICR96Pool = 1,
         Gene = Gene.Names,
         MLPAResult = "CNV",
         ResultType = "CNV",
         ExonCNVType = Type.of.CNV,
         ExonCNVSize = "CNV",
         Chromosome = Man..cur..Chr,
         `X5PrimeExon37` = Manual.curation.Start, 
         `X3PrimeExon37` = Manual.curation.Stop)

# --- 4. List of sample-gene pairs in TP rows ---
tp_samples_genes <- tp_rows %>% 
  select(sample, Gene.Names) %>% 
  mutate(SampleID = sample,
         Gene = Gene.Names) %>%
  distinct()

# --- 5. Generate full TN table and exclude TP rows ---
validated_results <- expand_grid(SampleID = sample_ids, dnascreen_genes) %>%
  rename(Gene = gene, Chromosome = chr, `X5PrimeExon37` = start, `X3PrimeExon37` = end) %>%
  anti_join(tp_samples_genes, by = c("SampleID", "Gene")) %>%
  mutate(
    ICR96Pool = 1,
    MLPAResult = "Normal",
    ResultType = "Normal",
    ExonCNVType = "",
    ExonCNVSize = ""
  ) %>%
  select(SampleID, ICR96Pool, Gene, MLPAResult, ResultType,
         ExonCNVType, ExonCNVSize, Chromosome, `X5PrimeExon37`, `X3PrimeExon37`)

# --- 6. Final result: TP + TNs ---
# Get the common column names
common_cols <- intersect(names(tp_rows), names(validated_results))

# Select only the common columns and bind them
final_validated <- bind_rows(
  tp_rows %>% select(all_of(common_cols)),
  validated_results %>% select(all_of(common_cols))
)

final_validated$SampleID <- paste0(final_validated$SampleID, ".hq.sorted.marked")

# --- 7. Save as TSV ---
# write.table(final_validated, file = "examples/final_validated_run_val_samples.tsv",
#             sep = "\t", row.names = FALSE, quote = FALSE)


