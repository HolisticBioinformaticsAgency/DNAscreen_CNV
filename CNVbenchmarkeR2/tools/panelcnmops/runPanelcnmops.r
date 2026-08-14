# Runs panelcn over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runpanelcn.R [panelcn_params_file] [datasets_params_file] [keepTempFiles]
# keepTempFiles: if true, temp files will not be removed (Default: true)
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions
suppressPackageStartupMessages(library(panelcn.mops))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(GenomicRanges))

#Functions----
# translates DEL/DUP into a common format
auxCNname <- function(x) {
  if (x %in% c("CN0", "CN1")) return("deletion")
  else if (x %in% c("CN3", "CN4")) return("duplication")
}

#Get parameters----
## Read args----
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
  panelcnParamsFile <- args[1]
  datasetsParamsFile <- args[2]
  keepTempFiles <- args[3]
} else {
  panelcnParamsFile <- "tools/panelcnmops/panelcnmopsParams.yaml"
  datasetsParamsFile <- "datasets.yaml"
  keepTempFiles <- "true"
}


##Load the parameters file----
params <- yaml.load_file(panelcnParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)

print(paste("Params for this execution:", list(params)))
print(paste("Datasets for this execution:", list(datasets)))

# Dataset iteration ----
# go over datasets and run panelcn for those which are active
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if (dataset$include){
    print(paste("Starting panelcn.mops for", name, "dataset", sep=" "))

    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bedFile <- file.path(dataset$bed_file)
    fastaFile <- file.path(dataset$fasta_file)
    validatedResultsFile <- file.path(dataset$validated_results_file)

    # set readlength from algorithm params if defined
    if (!is.null(params$readLength)){
      readLength <- params$readLength
    } else
      readLength <- dataset$read_length

    # Create output folder
    if (!is.null(params$outputFolder)) {
      outputFolder <- params$outputFolder
    } else
      outputFolder <- file.path(getwd(), "output", paste0("panelcnmops-", name))
    if (is.null(params$execution) || params$execution != "skipPrecalcPhase") {
      unlink(outputFolder, recursive = TRUE);
      dir.create(outputFolder)
    }
    outputFile <- file.path(outputFolder, "cnvFounds.txt")

    ##Do pre-calc part of the algorithm----
    if (is.null(params$execution) || params$execution != "skipPrecalcPhase") {
      # Get count windows
      countWindows <- getWindows(bedFile)

      # Get read counts from BAM file
      allbams <- list.files(path=bamsDir, pattern="*.bam$", full.names = TRUE)
      counts <- countBamListInGRanges(countWindows = countWindows, bam.files = allbams, read.width = readLength)
      names(mcols(counts)) <- unlist(lapply(names(mcols(counts)), getSampleName)) # remove .bam from samples names

      if (!is.null(params$execution) && params$execution == "onlyPrecalcPhase") {
        # save results to be used by other executions
        dir.create(params$precalcFolder)
        saveRDS(countWindows, file.path(params$precalcFolder, "countWindows.rds"))
        saveRDS(counts, file.path(params$precalcFolder, "counts.rds"))

        print(paste("panelcn.mops (Only pre-calc phase) for", name, "dataset finished", sep=" "))
        cat("\n\n\n")
        quit()
      }
    } else { # skipPrecalcPhase mode: read previous results
      print(paste("panelcn.mops Skipping pre-calc phase for", name, "dataset finished", sep=" "))
      countWindows <- readRDS(file.path(params$precalcFolder, "countWindows.rds"))
      counts <- readRDS(file.path(params$precalcFolder, "counts.rds"))
    }

    ##Extract params----
    classes <- c(params$CN0, params$CN1, 1, params$CN3, params$CN4)
    normType <- params$normType
    sizeFactor <- params$sizeFactor
    qu <- params$qu
    quSizeFactor <- params$quSizeFactor
    norm <- params$norm
    priorImpact <- params$priorImpact
    minMedianRC <- params$minMedianRC
    maxControls <- params$maxControls
    corrThresh  <- params$corrThresh

    ## Test samples depending on sample indications compatibility----
    allResults <- data.frame()
    allSampleNames <- names(mcols(counts))
    if (dataset$validated_results_file_format == "panelcn" && dataset$validated_results_file != ""
        && params$defineControlsByIndication == TRUE) {
      sample_indications <- readIndicationsForSamples(dataset$validated_results_file)
      indicationsComp <- readCompatibleIndications(dataset$validated_results_file) # get compatible indications for each indication

      for (indication in names(indicationsComp)) {
        print(paste("Processing indication:", indication))

        # samples to test: samples matching this indication
        testSamplesNames <- sample_indications[sample_indications$Genes == indication,][["SampleID"]]
        if (length(testSamplesNames > 0)){
          testSamples <- counts  # copy calculated counts
          mcols(testSamples) <- NULL # remove metadata columns
          mcols(testSamples)[, testSamplesNames] <- mcols(counts)[, testSamplesNames] # add only test samples

          # control samples: samples with compatible indications
          compatibles <- indicationsComp[[indication]]
          samplesToExclude <- sample_indications[!sample_indications$Genes %in% compatibles,][["SampleID"]]
          controlSamplesNames <- allSampleNames[!(allSampleNames %in% samplesToExclude)]
          controlSamples <- counts  # copy calculated counts
          mcols(controlSamples) <- NULL # remove metadata columns
          mcols(controlSamples)[, controlSamplesNames] <- mcols(counts)[, controlSamplesNames] # add only control samples

          ## Run the algorithm
          XandCB <- testSamples
          elementMetadata(XandCB) <- cbind(elementMetadata(XandCB), elementMetadata(controlSamples))
          resultList <- runPanelcnMops(XandCB, 1:ncol(elementMetadata(testSamples)),countWindows = countWindows, I = classes, sizeFactor = sizeFactor, norm = norm,
                                       normType = normType, qu = qu, quSizeFactor = quSizeFactor, priorImpact = priorImpact, minMedianRC = minMedianRC, maxControls = maxControls, corrThresh = corrThresh )

          # Build results table
          sNames <- colnames(elementMetadata(testSamples))
          resultsTable <- createResultTable(resultlist = resultList, XandCB = XandCB, countWindows = countWindows,
                                            sampleNames = sNames)
          concatResults <- ldply(resultsTable, data.frame) # concat output from all samples

          # join with final table
          allResults <- rbind(allResults, concatResults)
        }
      }
    } else {
      message("Clinical indication not found, using all samples potentially as controls")

      ## Run the algorithm----
      XandCB <- counts
      elementMetadata(XandCB) <- cbind(elementMetadata(XandCB), elementMetadata(XandCB))
      resultList <- runPanelcnMops(XandCB, 1:ncol(elementMetadata(counts)),countWindows = countWindows, I = classes, sizeFactor = sizeFactor, norm = norm,
                                   normType = normType, qu = qu, quSizeFactor = quSizeFactor, priorImpact = priorImpact, minMedianRC = minMedianRC, maxControls = maxControls, corrThresh = corrThresh)


      # Build results table
      sampleNames <- colnames(elementMetadata(counts))
      finalResultsTable <- createResultTable(resultlist = resultList, XandCB = XandCB, countWindows = countWindows,
                                             sampleNames = sampleNames)
      allResults <- ldply(finalResultsTable, data.frame) # concat output from all samples
    }

    #Save results----
    ## Build output file----
    colNames <- c("Sample", "Gene", "Chr", "Start", "End", "lowQual", "CN")
    filteredResults <- allResults[(allResults$CN != "CN2") & (allResults$lowQual != "lowQual"),colNames] # only deletions/duplications (CN2 means normal), and high quality
    filteredResults$CNV.type <- lapply(filteredResults$CN, function(x) sapply(x, auxCNname)) # Add CNV.type column before storing file
    filteredResults$CNV.type <- as.factor(unlist(filteredResults$CNV.type))  # R things...
    #Edit so that it matches validation results
    filteredResults$Sample <- sub("\\.bam$", "", filteredResults$Sample)
    write.table(filteredResults, outputFile, sep="\t", row.names=FALSE, quote = FALSE)  # write output file

    ## Save failed ROIs in a common format----
    failedROIs <- allResults[allResults$lowQual == "lowQual", colNames] # get all low qual
    names(failedROIs)[1]<- "SampleID" # rename Sample column
    failedROIs <- failedROIs[,c(1,3,4,5,2)] # reorder and filter columns
    failedROIs[,1] <- unlist(strsplit(data.frame(lapply(failedROIs, as.character), stringsAsFactors=FALSE)[,1],"\\.bam"))  # remove .bam from sample names
    write.table(failedROIs, file.path(outputFolder, "failedROIs.csv"), sep="\t", row.names=FALSE, quote = FALSE) # save

    #GenomicRanges object
    message("Saving GenomicRanges results")
    saveResultsFileToGR(outputFolder, "cnvFounds.txt")

    print(paste("panelcn.mops for", name, "dataset finished", sep=" "))
    cat("\n\n\n")

    #Temporary files----
    #Delete temporary files if specified
    if(keepTempFiles == "false"){
      filesAll <- list.files(outputFolder, full.names = TRUE)
      filesToKeep <- c("failedROIs.csv", "grPositives.rds", "cnvs_summary.tsv", "cnvFounds.csv", "cnvFounds.txt", "all_cnv_calls.txt", "calls_all.txt", "failures_Failures.txt", "cnv_calls.tsv")
      filesToRemove <- list(filesAll[!(filesAll %in% grep(paste(filesToKeep, collapse = "|"), filesAll, value = TRUE))])
      do.call(unlink, filesToRemove)
    }
  }
}

print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)
