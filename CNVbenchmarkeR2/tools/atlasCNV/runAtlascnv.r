# Runs Atlas-CNV over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runAtlasCNV.R [atlasCNVparams_file] [datasets_params_file] [keepTempFiles] [evaluateParameters]
# keepTempFiles: if true, temp files will not be removed (Default: true)
# evalulateParameters: if true, DoC will not be computed (Default: false)

print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(gtools)

#Get parameters----
## Read args----
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
  paramsFile <- args[1]
  datasetsParamsFile <- args[2]
  keepTempFiles <- args[3]
  evaluateParameters <- args[4]
} else {
  paramsFile <- "params.yaml"
  datasetsParamsFile <- "../../datasets.yaml"
  keepTempFiles <- "true"
  evaluateParameters <- "false"
}

## Load the parameters file----
params <- yaml.load_file(paramsFile)
datasets <- yaml.load_file(datasetsParamsFile)

# print params
print(paste("Params for this execution:", list(params)))
print(paste("Datasets for this execution:", list(datasets)))

## Get AtlasCNV folder ----
atlascnvFolder <- file.path(params$atlascnvFolder)

## Locate the reference ----
gatkFolder <- file.path(params$gatkFolder)
currentFolder <- getwd()

# Dataset iteration ----
# go over datasets and run AtlasCNV for those which are active
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if (dataset$include){
    print(paste("Starting Atlas-CNV for", name, "dataset", sep=" "))

    # Create output folder
    if (!is.null(params$outputFolder)) {
      if(stringr::str_detect(params$outputFolder, "^./")) params$outputFolder <- stringr::str_sub(params$outputFolder, 3, stringr::str_length(params$outputFolder))
      outputFolder <- file.path(currentFolder, params$outputFolder)
    } else {
      outputFolder <- file.path(getwd(), "output", paste0("atlasCNV-", name))
    }

    unlink(outputFolder, recursive = TRUE);
    dir.create(outputFolder, showWarnings = FALSE)

    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)
    bedFile <- file.path(dataset$bed_file)
    fastaFile <- file.path(dataset$fasta_file)

    ## Input files for Atlas-CNV ----

    ### Panel file ----
    #create panel file to input in Atlas
    tempBed<-read.csv(bedFile, sep="\t", header = FALSE)

    panel <- tempBed %>%
      mutate(V1 = factor(V1, levels = mixedsort(unique(V1)))) %>%  # Sort chromosomes naturally
      arrange(V1, as.numeric(V2)) %>%  # Sort by start position numerically
      mutate(Exon_Target = paste0(V1, ":", as.numeric(V2) + 1, "-", V3)) %>%
      mutate(Gene_Exon = V4) %>%
      mutate(Call_CNV = ifelse(str_starts(Exon_Target, "Y"), "N", "Y")) %>%
      mutate(RefSeq = V4) %>%
      select(-starts_with("V"))
    panel$Gene_Exon <- ave(panel$Gene_Exon, panel$Gene_Exon, FUN = function(i) paste0(i, '_', seq_along(i)))

    #export panel file
    write.table(panel, file = paste0(outputFolder,'/panel',name,'.txt'),sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

    ### Sample file ----

    tempSample <- as.data.frame(sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(bamFiles)))
    sample <- tempSample %>%
      dplyr::mutate(sex = "F") %>%
      dplyr::mutate(mp = name)

    #export sample file
    write.table(sample, file = paste0(outputFolder,'/', name, '.sample'), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    #Remove intermediate files
    rm(tempSample, tempBed)

    ### Depth of Coverage----
    #create Depth of coverage folder (DOC)
    depthCoverageFolder <- ifelse(evaluateParameters == "false",
                                  file.path(outputFolder, "DepthOfCoverage"),
                                  paste0(currentFolder, "/evaluate_parameters/atlasCNV/", name, "/DepthOfCoverage"))
    #Check if DOC folder exists to skip this step
    if(!dir.exists(depthCoverageFolder)| dir.exists(depthCoverageFolder) &&  length(list.files(depthCoverageFolder))<7){
    # #if it doesn't exist it calculates DOC
    dir.create(depthCoverageFolder, showWarnings = FALSE)

    #get depth of coverage for each BAM file, using GATK
    bam <- basename(bamFiles) %>% tools::file_path_sans_ext()
    for (i in seq_len(length(bamFiles))){
      cmd <- paste("singularity exec -B /fs03/vh83 -B /fs04/scratch2/vh83", gatkFolder, "gatk DepthOfCoverage",
               "-R", fastaFile,
               "-I", bamFiles[i],
               "-O", file.path(depthCoverageFolder, paste0(bam[i], ".DATA")),
               "--output-format", "TABLE",
               "-L", bedFile)
      print(cmd);system(cmd)
      }
    }

    ##Run atlasCNV----
    print("run atlasCNV")
    setwd(outputFolder)

    #create configuration file
    ##atlasCNV configuration file for Linux users
    cfg <- data.frame(V1 = paste0("GATKDIR=", depthCoverageFolder, "/[SAMPLE_FCLBC].DATA.sample_interval_summary")) %>%
      dplyr::add_row(V1 = paste0("ATLASCNV=", params$atlascnvFolder)) %>%
      dplyr::add_row(V1 = paste0("RPATH=", params$rpath)) %>%
      dplyr::add_row(V1 = paste0("RSCRIPT=", params$Rscript))

    write.table(cfg, paste0(outputFolder,"/config"), col.names = FALSE, row.names = FALSE, quote = FALSE)


    cmd_run_atlascnv <- paste0( "perl ",atlascnvFolder,"/atlas_cnv.pl",
                                " --config ", outputFolder,"/config",
                                " --panel ", outputFolder,'/panel',name,'.txt',
                                " --sample ", outputFolder,"/", name, ".sample",
                                " --threshold_del ", params$threshold_del,
                                " --threshold_dup ", params$threshold_dup)

    print(cmd_run_atlascnv);system(cmd_run_atlascnv);
    # Save results----
    ##TSV file----
    #merge cnv output files and modify columns
    resDF <- list.files(name, pattern = "\\.cnv", recursive = TRUE, full.names = TRUE)  %>%
      purrr::set_names(basename) %>%
      purrr::map(read.delim) %>%
      rlist::list.rbind() %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "sample") %>%
      dplyr::mutate(sample = stringr::str_replace_all(sample, ".cnv.FAILED_sampleQC_and_sampleANOVA.[0-9]+|.cnv.FAILED_sampleANOVA.[0-9]+|.cnv..FAILED_sampleQC.[0-9]+|.cnv.FAILED_sampleQC.[0-9]|.cnv.FAILED_sampleANOVA|.cnv.[0-9]+|.cnv","")) %>%
      dplyr::rename( "gene" = "Gene_Exon") %>%
      tidyr::separate(Exon_Target, c("chr", "start","end")) %>%
      dplyr::mutate(CNV.type = ifelse(cnv == "del",
                                      "deletion",
                                      "duplication"))

    # Read output file, add CNV.type column and write the results in a tsv file
    write.table(resDF, file.path(outputFolder, "cnv_calls.tsv"), sep="\t", quote=F, row.names = FALSE, col.names = TRUE)


    # Path to  tsv file
    finalSummaryFile <- file.path(outputFolder, "cnv_calls.tsv")
    ## GenomicRanges object ----
    message("Saving CNV GenomicRanges results")
    saveResultsFileToGR(outputFolder, basename(finalSummaryFile), geneColumn = "gene",
                        sampleColumn = "sample", chrColumn = "chr", startColumn = "start",
                        endColumn = "end", cnvTypeColumn = "CNV.type")

    ## Temporary files----
    #Delete temporary files if specified
    if(keepTempFiles == "false"){
      filesAll <- list.files(outputFolder, full.names = TRUE, recursive = TRUE)
      filesToKeep <- c("failedRois.csv", "grPositives.rds", "cnvs_summary.tsv", "cnvFounds.csv", "cnvFounds.txt", "all_cnv_calls.txt", "calls_all.txt", "failures_Failures.txt", "cnv_calls.tsv")
      filesToRemove <- list(filesAll[!(filesAll %in% grep(paste(filesToKeep, collapse = "|"), filesAll, value = TRUE))])
      do.call(unlink, filesToRemove)
    }
  }
#return to the original folder
  setwd("../..")
}

print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)

