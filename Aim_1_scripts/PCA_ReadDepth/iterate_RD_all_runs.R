library(ExomeDepth)
library(tools)
library(dplyr)
library(ggplot2)

# Record start time
start_time <- Sys.time()

with_Rdata = TRUE

dnascreen_runs = c(4:35)

# File paths
Rdata_dir <- "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/exome_depth/RData"
directory_bam <- paste0("/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/exome_depth/bams_run", dnascreen_run)
directory_vcf <- paste0("/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/VarSeq_CNV/vcfs_run", dnascreen_run)

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

# Filter bam_files based on existing sample IDs
filtered_bam_files <- bam_files[passqc_samples]

filtered_bam_filenames <- basename(filtered_bam_files)
filtered_bam_filenames <- gsub("-", ".", filtered_bam_filenames)

ref_fa = '/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/exome_depth/Homo_sapiens_assembly38.fasta'
# bed_name = "9genes_25bp.fix.sorted.bed"
bed_name = "3427651_Covered.bed"

bed_file = paste0('/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/CNV_calling/exome_depth/bed_files/', bed_name)


# Iterate through each run (4 to 35)
results_list <- list()
sample_ED_corr <- data.frame(
  sample = character(),
  run = character(),
  ref_test_corr = numeric(),
  stringsAsFactors = FALSE
)

for (run in dnascreen_runs) {
  cat("\nProcessing Run:", run, "\n")
  print(Rdata_dir)
  # Load read counts from Run
  
  if (with_Rdata){
    RData_file <- paste0(Rdata_dir, "/run", run,"_bams_9genes_25bp.fix.sorted.RData")
    load(RData_file)
  } else {
    my.count <- getBamCounts(bed.file = bed_file,
                             bam.files = filtered_bam_files,
                             include.chr = FALSE,
                             referenceFasta = ref_fa)
    save(my.count, file = paste0("RData/run", dnascreen_run, "_bams_", bed_name ,".RData"))
  }
 
  
  my.count.dafr <- as(my.count, 'data.frame')
  
  sample_ids <- colnames(my.count[6:length(my.count)])
  
  
  # Run CNV calling
  for (test_sample_id in sample_ids) {
    print(paste0('Running for ', test_sample_id))
    
    # Extract test sample counts
    my.test <- my.count[[test_sample_id]]
    
    # check if there are reads in sample
    if (sum(my.test == 0) > 10) {
      cat("Skipping", test_sample_id, "due to low coverage\n")
      next
    }
    
    # making sure the ref sample pool is without the test sample
    my.ref.samples <- sample_ids[!sample_ids %in% test_sample_id]
    my.reference.set <- as.matrix(my.count.dafr[, my.ref.samples])
    my.choice <- select.reference.set (test.counts = my.test,
                                       reference.counts = my.reference.set,
                                       bin.length = (my.count.dafr$end - my.count.dafr$start)/1000,
                                       n.bins.reduced = 10000)
    # Remove duplicates and additional suffix
    my.choice$reference.choice <- unique(gsub("\\.1$", "", my.choice$reference.choice))
    
    
    my.matrix <- as.matrix(my.count.dafr[, my.choice$reference.choice, drop = FALSE])
    my.reference.selected <- apply(X = my.matrix,
                                   MAR = 1,
                                   FUN = sum)
    
    all.panel <- new('ExomeDepth',
                     test = my.test,
                     reference = my.reference.selected,
                     formula = 'cbind(test, reference) ~ 1')
    
    output <- capture.output({
      all.panel <- CallCNVs(x = all.panel,
                            transition.probability = 10^-4,
                            chromosome = my.count.dafr$chromosome,
                            start = my.count.dafr$start,
                            end = my.count.dafr$end,
                            name = my.count.dafr$exon)
    }, type = "message")
    
    # Find the line containing the correlation
    cor_line <- grep("Correlation between reference and tests count", output, value = TRUE)
    
    # Extract the numeric correlation value
    correlation <- as.numeric(sub(".*is ([0-9\\.]+).*", "\\1", cor_line))
    
    # Store sample ref_test_corr in sample_ED_corr
    normal_sample_id <- gsub("\\.", "-", gsub("\\.hq\\.sorted\\.marked\\.bam$", "", test_sample_id))
    new_row <- data.frame(sample = normal_sample_id, run = run, ref_test_corr = correlation, stringsAsFactors = FALSE)
    sample_ED_corr <- rbind(sample_ED_corr, new_row)
    
    # Store results
    if (nrow(all.panel@CNV.calls) > 0) {
      all.panel@CNV.calls$Test_Filename <- test_sample_id
      results_list[[test_sample_id]] <- all.panel@CNV.calls
      results_list[[test_sample_id]]$Run <- run
      results_list[[test_sample_id]]$Correlation_ref_test_count <- correlation
    }
  }
}


# Combine results
if (length(results_list) > 0) {
  combined_df <- do.call(rbind, results_list)
  # combined_df$Test_Filename <- gsub("\\.", "-", combined_df$Test_Filename)
  combined_df$sample_id <- gsub("\\.", "-", combined_df$Test_Filename)
  combined_df$sample_id <- sub("(_S\\d+).*", "", combined_df$sample_id)
  
  # Save results
  output_file = "combined_calls_ED_all_runs_v3.csv"
  write.csv(combined_df, output_file, row.names = FALSE)
  
  cat(paste0("Results saved to", output_file, "\n"))
} else {
  cat("No CNVs detected for any samples.\n")
}

# Record and print execution time
end_time <- Sys.time()
cat("Total execution time:", end_time - start_time, "\n")

########################################################################################################################

## Use the ExDP function for the plots for the samples called for CNVs.

