# Runs clearCNV over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runclearCNV.R [clearCNV_params_file] [datasets_params_file] [keepTempFiles]
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions
library(dplyr)

#Get parameters----
## Read args ----
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
  clearCNVParamsFile <- args[1]
  datasetsParamsFile <- args[2]
  keepTempFiles <- args[3]
} else {
  clearCNVParamsFile <- "tools/clearCNV/clearCNVParams.yaml"
  datasetsParamsFile <- "datasets.yaml"
  keepTempFiles  <- "true"
}

##Load the parameters files----
params <- yaml.load_file(clearCNVParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)
print(paste("Params for this execution:", list(params)))

##Get clearCNV folder----
clearCNVFolder <- file.path(params$clearCNVFolder)

# Dataset iteration ----
# go over datasets and run clearCNV for those which are active
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if (dataset$include){
    print(paste("Starting clearCNV for", name, "dataset", sep=" "))

    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bedFile <- file.path(dataset$bed_file)
    fastaFile <- file.path(dataset$fasta_file)


    # Create output folder
    if (!is.null(params$outputFolder)) {
      outputFolder <- params$outputFolder
    } else {
      outputFolder <- file.path(getwd(), "output", paste0("clearCNV-", name))
    }
    unlink(outputFolder, recursive = TRUE);
    dir.create(outputFolder, showWarnings = FALSE)
    ## Input files for clearCNV ----
    ### Get bam files
    bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)
    ## Create bam directory and merge all bams in a txt file
    dir.create(file.path(outputFolder, "bams"))
    bamTxt <- file.path(outputFolder, paste0("bams/all_bams.txt"))
    write.table(x = paste(bamFiles, sep = "\n"), file = bamTxt , quote = FALSE, row.names = FALSE, col.names = FALSE)


    # Run clearCNV ----
    cmd <- paste(
        "source activate mamba-env &&",  # Activate mamba-env
        "clearCNV",
        "workflow_cnv_calling",
        "-w", outputFolder,
        "-p", name,
        "-r", fastaFile,
        "-b", bamTxt,
        "-d", bedFile,
        "-k", file.path("tools/clearCNV", "consensusBlacklist.bed"),
        "-c", params$cores,
        "--expected_artefacts", params$expected_artefacts,
        "--sample_score_factor", params$sample_score_factor,
        "--minimum_group_sizes", params$minimum_group_sizes,
        "--zscale", params$zscale,
        "--size", params$size,
        "--del_cutoff", params$del_cutoff,
        "--dup_cutoff", params$dup_cutoff,
        "--trans_prob", params$trans_prob
    )
    paste(cmd);system(cmd);

    # Read output file, add CNV.type column and write the results in a tsv file
    resFile <- file.path(outputFolder, paste0(name, "/results/cnv_calls.tsv"))
    resDF <- read.table(resFile, header = TRUE) %>%
      dplyr::mutate(CNV.type = ifelse(aberration == "DEL",
                                      "deletion",
                                      "duplication"))
    write.table(resDF, file.path(outputFolder, "cnv_calls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

    # Save results----
    # Path to  tsv file
    finalSummaryFile <- file.path(outputFolder, "cnv_calls.tsv")
    ## GenomicRanges object ----
    message("Saving CNV GenomicRanges results")
    saveResultsFileToGR(outputFolder, basename(finalSummaryFile), geneColumn = "gene",
                        sampleColumn = "sample", chrColumn = "chr", startColumn = "start",
                        endColumn = "end", cnvTypeColumn = "CNV.type")

    ## Temporary files----
    #Delete temporary files if specified
    if(keepTempFiles  == "false"){
      filesAll <- list.files(outputFolder, full.names = TRUE, recursive = TRUE)
      filesToKeep <- c("failedROIs.csv", "grPositives.rds", "cnvs_summary.tsv", "cnvFounds.csv", "cnvFounds.txt", "all_cnv_calls.txt", "calls_all.txt", "failures_Failures.txt", "cnv_calls.tsv")
      filesToRemove <- list(filesAll[!(filesAll %in% grep(paste(filesToKeep, collapse = "|"), filesAll, value = TRUE))])
      do.call(unlink, filesToRemove)
    }

  }
}
