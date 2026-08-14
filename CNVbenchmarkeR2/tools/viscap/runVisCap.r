# Runs viscap over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runVisCap.R [viscap_params_file] [datasets_params_file] [keepTempFiles] [evaluateParameters]
# keepTempFiles: if true, temp files will not be removed (Default: true)
# evalulateParameters: if true, DoC will not be computed (Default: false)
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions
library(dplyr)
library(stringr)


#Get parameters----
## Read args----
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
  viscapParamsFile <- args[1]
  datasetsParamsFile <- args[2]
  keepTempFiles <- args[3]
  evaluateParameters <- args[4]
} else {
  viscapParamsFile <- "tools/viscap/viscapParams.yaml"
  datasetsParamsFile <- "datasets.yaml"
  keepTempFiles <- "true"
  evaluateParameters <- "false"
}

# # remove this once you fixed R library issue and use original runBenchmark.R script
# evaluateParameters <- "false"


##Load the parameters file----
params <- yaml.load_file(viscapParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)
# print params
print(paste("Params for this execution:", list(params)))
print(paste("Datasets for this execution:", list(datasets)))


#Get viscap folder names----
currentFolder <- getwd()
viscapFolder <- file.path(params$viscapFolder)
gatkFolder <- file.path(params$gatkFolder)


# Dataset iteration ----
# go over datasets and run viscap for those which are active
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if (dataset$include){
    print(paste("Starting viscap for", name, "dataset", sep=" "))

    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bedFile <- file.path(dataset$bed_file)
    fastaFile <- file.path(dataset$fasta_file)


    # Create output folder
    if (!is.null(params$outputFolder)) {
        if(stringr::str_detect(params$outputFolder, "^./")) params$outputFolder <- stringr::str_sub(params$outputFolder, 3, stringr::str_length(params$outputFolder))
        outputFolder <- file.path(currentFolder, params$outputFolder)
    } else {
      outputFolder <- file.path(getwd(), "output", paste0("viscap-", name))
    }
    unlink(outputFolder, recursive = TRUE);
    dir.create(outputFolder, showWarnings = FALSE)



    #create input files required for running GATK

    ##Depth of coverage (DOC)----
    bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)
    depthCoverageFolder <- ifelse(evaluateParameters=="false",
                                  file.path(outputFolder, "DepthOfCoverage"),
                                  paste0(currentFolder, "/evaluate_parameters/viscap/", name, "/DepthOfCoverage")

    )
    #if DOC exists this step is skipped
    if(!dir.exists(depthCoverageFolder)| dir.exists(depthCoverageFolder) &&  length(list.files(depthCoverageFolder))<7){
    dir.create(depthCoverageFolder, showWarnings = FALSE)
    bam <- basename(bamFiles) %>% tools::file_path_sans_ext()
    #Create DOC files for all samples
    for (i in seq_len(length(bamFiles))){
      cmd <- paste("singularity exec -B /fs03/vh83 -B", paste0(dirname(fastaFile), ",", outputFolder),  
        gatkFolder, "gatk DepthOfCoverage",
                   "-R",  fastaFile,
                   "-I", bamFiles[i],
                   "-O",  file.path(depthCoverageFolder,bam[i]),
                   "--output-format", "TABLE",
                   "-L",  bedFile )
      paste(cmd);system(cmd)
    }
    }
    #need to enter to the folder where viscap config is stored
    setwd(outputFolder)

    ##Create a cfg file with all the parameters to use----
    viscapConfig<- file.path(outputFolder, "VisCap.cfg")
    cfg <- data.frame(V1=c("#VisCap configuration file for Linux users"))
    cfg <- cfg %>%
      dplyr::add_row(V1=paste0('interval_list_dir           <- "', dirname(bedFile), '"' )) %>%
      dplyr::add_row(V1=paste0('cov_file_pattern            <- "', params$cov_file_pattern, '"'  )) %>%
      dplyr::add_row(V1=paste0('cov_field                   <- "', params$cov_field, '"' )) %>%
      dplyr::add_row(V1=paste0('interval_file_pattern       <- "^', basename(bedFile), '$"' )) %>%
      dplyr::add_row(V1=paste0('ylimits            <- ', params$ylimits)) %>%
      dplyr::add_row(V1=paste0('iqr_multiplier              <- ', params$iqr_multiplier)) %>%
      dplyr::add_row(V1=paste0('threshold.min_exons  <- ', params$threshold.min_exons)) %>%
      dplyr::add_row(V1=paste0('threshold.cnv_log2_cutoffs  <- ', 'c(', params$threshold.cnv_log2_cutoffs_min, ', ', params$threshold.cnv_log2_cutoffs_max, ')' )) %>%
      dplyr::add_row(V1=paste0('iterative.calling.limit     <- ', params$iterative.calling.limit)) %>%
      dplyr::add_row(V1=paste0('infer.batch.for.sub.out_dir  <- ', params$infer.batch.for.sub.out_dir)) %>%
      dplyr::add_row(V1=paste0('clobber.output.directory    <- ', params$clobber.output.directory)) %>%
      dplyr::add_row(V1=paste0('dev_dir     <- ', params$dev_dir))
    write.table(cfg,  viscapConfig, col.names = FALSE, row.names = FALSE, quote=FALSE)

    if(params$iterative.calling.limit==1){
      outputFolder1 <- file.path(outputFolder, paste0(basename(outputFolder), "_run1"))
      dir.create(outputFolder1)
    }else{
      outputFolder1 <- outputFolder
    }

    ##Run Viscap----
    cmd <- paste("Rscript", file.path(viscapFolder, "VisCap.R"),
                 depthCoverageFolder,
                 outputFolder1,
                 viscapConfig)

    print(cmd); system(cmd)

    setwd("../..") #return to the original folder

    #Create a dataframe to store viscap results
    cnvData <- data.frame()
    ## Get the result of the most recent run
    runs <- list.files(path=outputFolder, pattern="_run[0-9]+$", full.names=TRUE)
    numberRuns <- length(runs)
    run <- runs[numberRuns]
    resultsFiles <- list.files(path=run, pattern=".cnvs.xls", full.names = TRUE)
    for(i in seq_len(length(resultsFiles))){
      testResultDf <- read.csv2(resultsFiles[i], sep="\t")%>%
        dplyr::mutate(chr= stringr::str_split(Genome_start_interval, ":") %>% purrr::map(1) %>% unlist(),
                      start = stringr::str_split(Genome_start_interval, ":|-") %>% purrr::map(2) %>% unlist(),
                      end= stringr::str_split(Genome_end_interval, ":|-") %>% purrr::map(3) %>% unlist()) %>%
        dplyr::mutate(copy_number = ifelse(CNV == "Loss",
                                           "deletion",
                                           ifelse(CNV=="Gain",
                                                  "duplication",
                                                  "")))

      cnvData <- rbind(cnvData, testResultDf)
    }
    cnvData <- cnvData %>% select("chr", "start", "end", "copy_number", "Median_log2ratio", "Interval_count", "Sample" )
    names(cnvData) <- c("chr", "start", "end", "copy_number", "quality", "targets", "Sample")

    #Convert bed file and cnvData results to GRanges class object
    bedGR <- regioneR::toGRanges(bedFile)
    resultsGR <- regioneR::toGRanges(cnvData)
    # Create a dataframe to store results with the associated gene
    resGenes <- data.frame()
    # go over each CNV to identify genes
    for (i in seq_len(length(resultsGR))){

      #get the CNV
      cnvGR <- regioneR::toGRanges(resultsGR[i,])

      # Merge the cnv with the bed file to know the "affected" genes.
      # Mutate columns to include information related to the cnv of interest (quality, the copy_number, Sample..)
      #Group by gene to only have one row per gene and get the smallest start coordinate and the biggest end coordinate
      res <- as.data.frame(IRanges::subsetByOverlaps(bedGR, cnvGR)) %>%
        dplyr::mutate(copy_number = cnvGR$copy_number,
                      quality = cnvGR$quality,
                      targets = cnvGR$targets,
                      Sample = resultsGR$Sample[i]) %>%
        dplyr::group_by(name, Sample, seqnames) %>%
        dplyr::summarise(start_gene = min(start),
                         end_gene= max (end),
                         copy_number = min(copy_number),
                         quality = min(quality),
                         width = sum(width))

      # Merge each cnv result per gene with the rest
      resGenes <- rbind(resGenes, res)
    }

    # Save results----
    ## TSV file----
    finalSummaryFile <- file.path(outputFolder, "cnvs_summary.tsv")
    write.table(resGenes, finalSummaryFile, sep = "\t", quote = F, row.names = F)

    ## GenomicRanges object----
    message("Saving CNV GenomicRanges results")
    saveResultsFileToGR(outputFolder, basename(finalSummaryFile), geneColumn = "name",
                        sampleColumn = "Sample", chrColumn = "seqnames", startColumn = "start_gene",
                        endColumn = "end_gene", cnvTypeColumn = "copy_number")

    ## Temporary files----
    #Delete temporary files if specified
    if(keepTempFiles == "false"){
      filesAll <- list.files(outputFolder, full.names = TRUE, recursive=TRUE)
      filesToKeep <- c("failedROIs.csv", "grPositives.rds", "cnvs_summary.tsv", "cnvFounds.csv", "cnvFounds.txt", "all_cnv_calls.txt", "calls_all.txt", "failures_Failures.txt", "cnv_calls.tsv")
      filesToRemove <- list(filesAll[!(filesAll %in% grep(paste(filesToKeep, collapse= "|"), filesAll, value=TRUE))])
      do.call(unlink, filesToRemove)
    }
  }
}