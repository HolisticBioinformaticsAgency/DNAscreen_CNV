# Runs cnvkit over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runCnvkit.R [cnvkit_params_file] [datasets_params_file] [keepTempFiles]
# keepTempFiles: if true, temp files will not be removed (Default: true)
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions
library(dplyr)

#Get parameters----
## Read args ----
args <- commandArgs(TRUE)
if(length(args)>0) {
  cnvkitParamsFile <- args[1]
  datasetsParamsFile <- args[2]
  keepTempFiles <- args[3]
} else {
  cnvkitParamsFile <- "tools/cnvkit/cnvkitParams.yaml"
  datasetsParamsFile <- "datasets.yaml"
  keepTempFiles <- "true"
}

##Load the parameters files----
params <- yaml.load_file(cnvkitParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)

# print params
print(paste("Params for this execution:", list(params)))
print(paste("Datasets for this execution:", list(datasets)))

##Get cnvkit folder----
cnvkitFolder <- file.path(params$cnvkitFolder)

# Dataset iteration ----
# go over datasets and run cnvkit for those which are active
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if (dataset$include){
    print(paste("Starting cnvkit for", name, "dataset", sep=" "))

    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bedFile <- file.path(dataset$bed_file)
    fastaFile <- file.path(dataset$fasta_file)
    currentFolder <- getwd()

    # Create output folder
    if (!is.null(params$outputFolder)) {
      if(stringr::str_detect(params$outputFolder, "^./")) params$outputFolder <- stringr::str_sub(params$outputFolder, 3, stringr::str_length(params$outputFolder))
      outputFolder <- file.path(currentFolder, params$outputFolder)
    } else {
      outputFolder <- file.path(getwd(), "output", paste0("cnvkit-", name))
    }
    unlink(outputFolder, recursive = TRUE);
    dir.create(outputFolder, showWarnings = FALSE)

    # Get bam files
    bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)
    ## Run cnvkit----
    ### Access----
    #calculate accessible coordinates of the chromosome sequences
    cmd <- paste(".",
                 cnvkitFolder, "\n",
                 "cnvkit.py",  "access",
                 fastaFile,
                 "-s", params$minGapSizeAccess,
                 "-o", file.path(outputFolder, "access.bed"))

    print(cmd);system(cmd)

    ###Autobin----
    #Quickly estimate read counts or depths in a BAM file to estimate reasonable on- and (if relevant) off-target bin sizes.
    cmd <- paste(".",
                 cnvkitFolder, "\n",
                 "cnvkit.py", "autobin",
                 paste0(bamsDir, "/*.bam"),
                 "-m hybrid",
                 "-t", bedFile,
                 "-g", file.path(outputFolder, "access.bed"),
                 "-b", params$bpPerBinAutobin,
                 "--target-max-size", params$targetMaxSizeAutobin,
                 "--target-min-size", params$targetMinSizeAutobin,
                 "--antitarget-max-size", params$atargMaxSizeAutobin,
                 "--antitarget-min-size", params$atargMinSizeAutobin,
                 "--target-output-bed",  file.path(outputFolder, "ROIs-ICR96-panelcnDataset.target.bed"),
                 "--antitarget-output-bed",  file.path(outputFolder, "ROIs-ICR96-panelcnDataset.antitarget.bed")
    )
    print(cmd);system(cmd)
    #create a folder to store calls files
    dir.create(file.path(outputFolder, "calls"))


    for(i in seq_len(length(bamFiles))){
      #For each sample
      sampleName <- basename(tools::file_path_sans_ext(bamFiles))

      #Get the coverage for target and antitarget
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   "cnvkit.py",  "coverage",
                   bamFiles[i],
                   "-q", params$minMapqCoverage,
                   file.path(outputFolder, "ROIs-ICR96-panelcnDataset.target.bed"),
                   "-o", file.path(outputFolder, paste0(sampleName[i], ".targetcoverage.cnn" )))
      print(cmd);system(cmd)

      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   #"srun",
                   "cnvkit.py","coverage",
                   bamFiles[i],
                   "-q", params$minMapqCoverage,
                   file.path(outputFolder, "ROIs-ICR96-panelcnDataset.antitarget.bed"),
                   "-o", file.path(outputFolder, paste0(sampleName[i], ".antitargetcoverage.cnn" )))
      print(cmd); system(cmd)
      }
    ### Coverage of target and antitarget----
    for(i in seq_len(length(bamFiles))){
      #Build a reference for each sample, including all samples except the sample of interest
      #Get cnn Target and antitarget files
      cnnTarget <- paste0(file.path(outputFolder, paste0(sampleName[-i], ".targetcoverage.cnn" )), collapse = " ")
      cnnAntitarget <- paste0(file.path(outputFolder, paste0(sampleName[-i], ".antitargetcoverage.cnn" )), collapse=" ")
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   "cnvkit.py", "reference",
                   cnnTarget, cnnAntitarget,
                   "-f", fastaFile,
                   "--min-cluster-size", params$minClusterSizeReference,
                   "-o", file.path(outputFolder, paste0(sampleName[i], "reference.cnn" )))
      paste(cmd);system(cmd)



      ###Fix----
      #uncorrected target and antitarget coverage tables (.cnn) and correct for biases in regional coverage and GC content, according to the given reference.
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   "cnvkit.py", "fix",
                   file.path(outputFolder, paste0(sampleName[i], ".targetcoverage.cnn" )),
                   file.path(outputFolder, paste0(sampleName[i], ".antitargetcoverage.cnn" )),
                   file.path(outputFolder, paste0(sampleName[i], "reference.cnn" )),
                   "-o",  file.path(outputFolder, paste0(sampleName[i], ".cnr")))
      paste(cmd); system(cmd)


      ###Segment----
      #infer discrete copy number segments from the given coverage table:
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   "cnvkit.py", "segment",
                   file.path(outputFolder, paste0(sampleName[i], ".cnr")),
                   "-m", params$methodSegment,
                   #"-t", params$thresehold,
                   "--drop-outliers", params$DropOutliersSegment,
                   "-o", file.path(outputFolder, paste0(sampleName[i], ".cns")))
      paste(cmd);system(cmd)

      ###Segmetrics----
      #Calculate summary statistics of the residual bin-level log2 ratio estimates from the segment means, similar to the existing metrics command, but for each segment individually.
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   "cnvkit.py", "segmetrics",
                   file.path(outputFolder, paste0(sampleName[i], ".cnr")),
                   "-s", file.path(outputFolder, paste0(sampleName[i], ".cns")),
                   "--ci",
                   "--alpha", params$alphaSegmetrics,
                   "-b", params$bootstrapSegmetrics,
                   "-o", file.path(outputFolder, paste0(sampleName[i],"segment.cns")))

      paste(cmd);system(cmd)

      ###Call cnvs----
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   "cnvkit.py", "call",
                   file.path(outputFolder, paste0(sampleName[i], "segment.cns")),
                   "--filter", params$filterCall,
                   "-m", params$methodCall,
                   paste0("-t=", params$thresholdsCall_del, ",", params$thresholdsCall_loss, ",", params$thresholdsCall_gain, ",", params$thresholdsCall_amp),
                   "-o", file.path(outputFolder, paste0("calls/", sampleName[i], ".call.cns")))
      paste(cmd);system(cmd)

      ###Bintest----
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   "cnvkit.py", "bintest",
                   file.path(outputFolder, paste0(sampleName[i], ".cnr")),
                   "-s", file.path(outputFolder, paste0("calls/", sampleName[i], ".call.cns")),
                   "-t",
                   "-a", params$alphaBintest,
                   "-o", file.path(outputFolder, paste0(sampleName[i], "bintest.cns"))
      )
      print(cmd); system(cmd)

      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   "cnvkit.py", "call",
                   file.path(outputFolder, paste0(sampleName[i], "bintest.cns")),
                   "-m", params$methodCall,
                   #"-t", params$thresholdsCall,
                   "-o", file.path(outputFolder, paste0("calls/", sampleName[i], "bintest.call.cns")))

      print(cmd); system(cmd)
      }


    #Save results----
    ##TSV file----
    #Put all results in the same file
    cnvData  <- data.frame()
    for(s in sampleName){
      #read the call files and merge them
      # Read call files
      callFiles <- read.table(file.path(outputFolder, paste0("calls/", s, ".call.cns")), header = TRUE)
      file_path <- file.path(outputFolder, paste0("calls/", s, "bintest.call.cns"))
      
      if (file.exists(file_path)) {
        callBinFiles <- read.table(file_path, header = TRUE)
      } else {
        callBinFiles <- data.frame()  # Empty dataframe
      }

      # Add missing columns if necessary
      if (nrow(callFiles) > 0) {
        callFiles <- callFiles %>%
          dplyr::mutate(p_bintest = NA) %>%
          dplyr::relocate(p_bintest, .after = depth)
      }

      if (nrow(callBinFiles) > 0) {
        callBinFiles$probes <- NA
      }

      # Check if both dataframes are empty before proceeding
      if (nrow(callFiles) > 0 & nrow(callBinFiles) > 0) {
        sam <- rbind(callFiles, callBinFiles)
      } else if (nrow(callFiles) > 0) {
        sam <- callFiles
      } else if (nrow(callBinFiles) > 0) {
        sam <- callBinFiles
      } else {
        next  # Skip to the next iteration if both are empty
      }

      #Put sample names and the consequence of the cnv
      sam2 <- sam %>%
        dplyr::rowwise() %>%
        dplyr::filter(cn != 2) %>%
        dplyr::mutate(Sample= s,
                      CNV.type = ifelse(cn < 2 ,
                                    "deletion",
                                    ifelse(cn>2,
                                           "duplication",
                                            "")))
        #Merge all in the same file
        cnvData <- rbind(cnvData, sam2)
      }

    #Read the bed file and pass cnv results to a GRanges object
    bedGR <- regioneR::toGRanges(bedFile)
    resultsGR <- regioneR::toGRanges(cnvData %>%
                                       as.data.frame())

    # Create a dataframe to store results with the associated gene
    resGenes <- data.frame()
    for (i in seq_len(length(resultsGR))){
      #get the CNV
      cnvGR <- regioneR::toGRanges(resultsGR[i,])
      # Merge the cnv with the bed file to know the "affected" genes.
      # Mutate columns to include information related to the cnv of interest (quality, the copy_number, Sample..)
      # Group by gene to only have one row per gene and get the smallest start coordinate and the biggest end coordinate
      res <- as.data.frame(IRanges::subsetByOverlaps(bedGR, cnvGR)) %>%
        dplyr::mutate(CNV.type = cnvGR$CNV.type,
                      log2 = cnvGR$log2,
                      cn = cnvGR$cn,
                      Sample = resultsGR$Sample[i]) %>%
        dplyr::group_by(name, Sample, seqnames, CNV.type) %>%
        dplyr::summarise(start_gene = min(start),
                         end_gene = max (end),
                         cn = min(cn),
                         log2 = min(log2))
      # Merge each cnv result per gene with the rest
      resGenes <- rbind(resGenes, res)
      }

    # Print results in a tsv file
    finalSummaryFile <- file.path(outputFolder, "cnvs_summary.tsv")
    write.table(resGenes, finalSummaryFile, sep = "\t", quote = FALSE, row.names = FALSE)

    ##GenomicRanges object----
    message("Saving CNV GenomicRanges results")
    saveResultsFileToGR(outputFolder, basename(finalSummaryFile), geneColumn = "name",
                        sampleColumn = "Sample", chrColumn = "seqnames", startColumn = "start_gene",
                        endColumn = "end_gene", cnvTypeColumn = "CNV.type")

    ##Temporary files----
    #Delete temporary files if specified
    if(keepTempFiles == "false"){
      filesAll <- list.files(outputFolder, full.names = TRUE)
      filesToKeep <- c("failedROIs.csv", "grPositives.rds", "cnvs_summary.tsv", "cnvFounds.csv", "cnvFounds.txt", "all_cnv_calls.txt", "calls_all.txt", "failures_Failures.txt", "cnv_calls.tsv")
      filesToRemove <- list(filesAll[!(filesAll %in% grep(paste(filesToKeep, collapse = "|"), filesAll, value=TRUE))])
      do.call(unlink, filesToRemove)
    }
  }
}

print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)

