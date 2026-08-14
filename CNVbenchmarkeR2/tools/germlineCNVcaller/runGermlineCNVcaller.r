# Runs GermlineCNVcaller on the datasets configured at [datasets_params_file]
# USAGE: Rscript runGermlineCNVcaller.R [GermllineCNVcaller_params_file] [datasets_params_file] [include_temp_files] [skipPrecalcPhase]
# keepTempFiles: if true, temp files will not be removed (Default: true)
# skipPrecalcPhase: If true, Initial steps will be omitted until DetermineGermlineContigPloidy. Used when calling evaluate parameters. (Default: false)

print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
source("utils/utils.r") # Load utils functions
source("tools/germlineCNVcaller/germlineCNVcallerUtils.r") # Load GATK wrap functions
library(dplyr)

#Get parameters----
## Read args----
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
  paramsFile <- args[1]
  datasetsParamsFile <- args[2]
  keepTempFiles <- args[3]
  skipPrecalcPhase <- args[4]
} else {
  paramsFile <- "params.yaml"
  datasetsParamsFile <- "../../datasets.yaml"
  keepTempFiles <- "true"
  skipPrecalcPhase <- "false"
}

print(paste("skipPrecalcPhase: ", skipPrecalcPhase))


## Load the parameters file----
params <- yaml.load_file(paramsFile)
datasets <- yaml.load_file(datasetsParamsFile)

# print params
print(paste("Params for this execution:", list(params)))
print(paste("Datasets for this execution:", list(datasets)))

## Load paths----
picard <- params$picard
gatk <- params$gatk
currentFolder <- getwd()
contigFile <- file.path(currentFolder, "tools/germlineCNVcaller/contig-ploidy-prior.tsv")


# Dataset iteration ----
# go over datasets and run germlineCNVcaller
for (name in names(datasets)) {
  dataset <- datasets[[name]]

  if (dataset$include){
    print(paste("Starting GermlineCNVcaller for", name, "dataset", sep=" "))

    # check if outputfolder exists
    if (!is.null(params$outputFolder)) {
      if(stringr::str_detect(params$outputFolder, "^./")) params$outputFolder <- stringr::str_sub(params$outputFolder, 3, stringr::str_length(params$outputFolder))
      outputFolder <- file.path(currentFolder, params$outputFolder)
    } else {
      outputFolder <- file.path(getwd(), "output", paste0("germlineCNVcaller-", name))
    }
    unlink(outputFolder, recursive = TRUE);
    dir.create(outputFolder)

    # extract dataset fields
    bamsDir <- file.path(dataset$bams_dir)
    bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)
    bedFile <- file.path(dataset$bed_file)
    fastaFile <- file.path(dataset$fasta_file)
    fastaDict <- file.path(dataset$fasta_dict)


    ## Create Interval list ----
    interval_list <- ifelse(skipPrecalcPhase == "false",
                            paste0(outputFolder,"/list.interval_list"),
                            paste0(currentFolder, "/evaluate_parameters/germlineCNVcaller/", name, "/preCalculated/IntervalList/list.interval_list"))
    if(skipPrecalcPhase == "false"){
      createIntervalList(picard, bedFile, interval_list, fastaDict)
    }


    ## GATK: Preprocess intervals ----
    preprocessInterval_list <- ifelse(skipPrecalcPhase=="false",
                                      paste0(outputFolder,"/preprocessed_intervals.interval_list"),
                                      paste0(currentFolder, "/evaluate_parameters/germlineCNVcaller/", name, "/preCalculated/IntervalList/preprocessed_intervals.interval_list"))
    if(skipPrecalcPhase == "false"){
      preprocessIntervals(gatk, fastaFile, outputFolder, interval_list, preprocessInterval_list)
    }

    ## GATK: Annotate intervals----
    annotatedIntervals_list <- ifelse(skipPrecalcPhase == "false",
                                      paste0(outputFolder,"/annotated_intervals.interval_list"),
                                      paste0(currentFolder, "/evaluate_parameters/germlineCNVcaller/", name, "/preCalculated/IntervalList/annotated_intervals.interval_list"))
    if(skipPrecalcPhase == "false"){
      annotateIntervals(gatk, fastaFile, outputFolder, preprocessInterval_list, annotatedIntervals_list)
    }


    ## GATK: Collect read counts ----
    bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)
    readCountsFolder <- ifelse(skipPrecalcPhase == "false",
                                  file.path(outputFolder, "ReadCounts"),
                                  paste0(currentFolder, "/evaluate_parameters/germlineCNVcaller/", name, "/preCalculated/ReadCounts"))
    if(skipPrecalcPhase == "false"){
      if(!dir.exists(readCountsFolder) | dir.exists(readCountsFolder) && length(list.files(readCountsFolder)) != length(bamFiles)){
        collectReadCounts(gatk, bamFiles, outputFolder, preprocessInterval_list, readCountsFolder)
      }
    }


    #merge all hdf5 files for running DetermineGermlineContigPloidy
    hdf5s <- list.files(readCountsFolder, pattern = '*.hdf5$', full.names = TRUE)
    all_hdf5_names<- paste(hdf5s, collapse=" -I ")

    ## GATK: Filter intervals----
    filteredIntervals_list <- ifelse(skipPrecalcPhase == "false",
                                     paste0(outputFolder,"/filtered_intervals.interval_list"),
                                     paste0(currentFolder, "/evaluate_parameters/germlineCNVcaller/", name, "/preCalculated/IntervalList/filtered_intervals.interval_list"))
    if(skipPrecalcPhase == "false"){
      filterIntervals(gatk, outputFolder, preprocessInterval_list, all_hdf5_names, annotatedIntervals_list, filteredIntervals_list)
    }

    ##GATK: Determine contig ploidy ----
    contigPloidyDir <- ifelse(skipPrecalcPhase == "false",
                               file.path(outputFolder, "contigPloidy"),
                               paste0(currentFolder, "/evaluate_parameters/germlineCNVcaller/", name, "/preCalculated/contigPloidy"))
    if(skipPrecalcPhase == "false"){
      determineContigPloidy(gatk, contigPloidyDir, bamsDir, outputFolder, contigFile, all_hdf5_names, filteredIntervals_list)
    }

    ##GATK: Run GermlineCNVcaller ----
    cnvRaw <- file.path(paste0(outputFolder, "/cnvRaw"))
    if (!file.exists(cnvRaw)){
      dir.create(cnvRaw)
    } else {
      unlink(cnvRaw, recursive = TRUE)
      dir.create(cnvRaw, recursive = TRUE)
    }

    print(paste("Starting at", Sys.time(), "GermlineCNVCaller"))
    cmd_GermlineCNVCaller <-  paste("source activate mamba-env &&", # Activate mamba-env
                                    "singularity exec -B /fs03/vh83 -B /fs04/vh83 -B /fs04/scratch2/vh83 -B", paste0(bamsDir, ",", outputFolder),
                                    gatk, 'gatk GermlineCNVCaller',
                                     "--run-mode COHORT",
                                     "-L ", filteredIntervals_list,
                                     "--interval-merging-rule OVERLAPPING_ONLY",
                                     "--annotated-intervals ", annotatedIntervals_list,
                                     "--contig-ploidy-calls ", paste0(contigPloidyDir, "/contig-calls"),
                                     "-I ", all_hdf5_names,
                                     "-O ", cnvRaw,
                                     "--output-prefix cohort_run",
                                     #optional parameters
                                      "--p-active", params$pActive,
                                      "--p-alt", params$pAlt,
                                      "--sample-psi-scale", params$samplePsiScale,
                                      "--mapping-error-rate", params$mappingErrorRate,
                                      "--class-coherence-length", params$classCoherenceLength,
                                      "--cnv-coherence-length", params$cnvCoherenceLength,
                                      "--interval-psi-scale", params$intervalPsiScale,
                                      "--active-class-padding-hybrid-mode", params$activeClassPaddingHybridMode,
                                      "--adamax-beta-1", params$adamaxBeta1,
                                      "--adamax-beta-2", params$adamaxBeta2,
                                      "--caller-external-admixing-rate", params$callerExternalAdmixingRate,
                                      "--caller-internal-admixing-rate", params$callerInternalAdmixingRate,
                                      "--caller-update-convergence-threshold", params$callerUpdateConvergenceThreshold,
                                      "--convergence-snr-averaging-window", params$convergenceSnrAveragingWindow,
                                      "--convergence-snr-countdown-window", params$convergenceSnrCountdownWindow,
                                      "--convergence-snr-trigger-threshold", params$convergenceSnrTriggerThreshold,
                                      "--depth-correction-tau", params$depthCorrectionTau,
                                      "--enable-bias-factors", params$enableBiasFactors,
                                      "--init-ard-rel-unexplained-variance", params$initArdRelUnexplainedVariance,
                                      "--learning-rate", params$learningRate,
                                      "--log-emission-samples-per-round", params$logEmissionSamplesPerRound,
                                      "--log-emission-sampling-median-rel-error", params$logEmissionSamplingMedianRelError,
                                      "--log-emission-sampling-rounds", params$logEmissionSamplingRounds,
                                      "--log-mean-bias-standard-deviation", params$logMeanBiasStandardDeviation,
                                      "--num-samples-copy-ratio-approx", params$numSamplesCopyRatioApprox,
                                      "--num-thermal-advi-iters", params$numThermalAdviIters,
                                      "--copy-number-posterior-expectation-mode", params$copyNumberPosteriorExpectationMode,
                                      "--gc-curve-standard-deviation", params$gcCurveStandardDeviation,
                                      "--max-bias-factors", params$maxBiasFactors,
                                      "--max-calling-iters", params$maxCallingIters,
                                      "--max-advi-iter-first-epoch", params$maxAdviIterFirstEpoch,
                                      "--max-advi-iter-subsequent-epochs", params$maxAdviIterSubsequentEpochs,
                                      "--max-training-epochs", params$maxTrainingEpochs,
                                      "--min-training-epochs", params$minTrainingEpochs,
                                      "--disable-annealing", params$disableAnnealing)

    print(cmd_GermlineCNVCaller);system(cmd_GermlineCNVCaller);


    ##GATK: Run PostprocessingGermlineCNVcaller----
    cnvProcessed <- file.path(paste0(outputFolder, "/cnvProcessed"))
    if (!file.exists(cnvProcessed)){
      dir.create(cnvProcessed)
    } else {
      unlink(cnvProcessed, recursive = TRUE)
      dir.create(cnvProcessed, recursive = TRUE)
    }

    callsRaw <- list.dirs(path=paste0(cnvRaw,"/cohort_run-calls"),  full.names = TRUE, recursive= FALSE)
    ploidyRaw <- list.dirs(path=paste0(contigPloidyDir,"/contig-calls"),  full.names = TRUE, recursive= FALSE)

    print(paste("Starting at", Sys.time(), "Germline postProcess"))
    for (i in 1:length(callsRaw)){
      cmd_GermlineCNVCaller <- paste("source activate mamba-env &&", # Activate mamba-env
                                     "singularity exec -B /fs03/vh83 -B /fs04/vh83 -B /fs04/scratch2/vh83 -B", paste0(bamsDir, ",", outputFolder),
                                     gatk, 'gatk PostprocessGermlineCNVCalls',
                                     "--sample-index",  i-1,
                                     "--calls-shard-path", paste0(cnvRaw, "/cohort_run-calls") ,
                                     "--contig-ploidy-calls", paste0(contigPloidyDir, "/contig-calls"),
                                     "--model-shard-path", paste0(cnvRaw, "/cohort_run-model"),
                                     "--output-genotyped-intervals", paste0(cnvProcessed, "/SAMPLE_", i-1 ,"_genotyped_intervals.vcf"),
                                     "--output-genotyped-segments", paste0(cnvProcessed, "/SAMPLE_", i-1  ,"_genotyped_segments.vcf"),
                                     "--output-denoised-copy-ratios", paste0(cnvProcessed, "/SAMPLE_", i-1 ,"_denoised_copy_ratios.tsv"))
       print(cmd_GermlineCNVCaller);system(cmd_GermlineCNVCaller)
    }


    # Save results ----
    ##TXT file----
    # Put all results in the same file
    resultsFiles <- list.files(path=cnvProcessed, pattern="_genotyped_segments.vcf$", full.names = TRUE, recursive= TRUE)
    cnvData  <- data.frame()
    for(i in seq_len(length(resultsFiles))){
     vcf <-vcfR::read.vcfR(resultsFiles[i])
     sample.name <- (vcf@gt %>% colnames())[2]
     cnvs_proc <- vcfR::vcfR2tidy(vcf, single_frame = TRUE)$dat %>%
       dplyr::rowwise() %>%
       mutate(sample=sample.name)%>%
       dplyr::filter(ALT != ".") %>%
       dplyr::mutate(copy_number = ifelse(ALT == "<DEL>",  "deletion", ifelse(ALT=="<DUP>","duplication", "")),
                     start = POS %>% as.integer(), end = END %>% as.integer())
     cnvData <- rbind(cnvData, cnvs_proc)
    }
    # so that it matches with the validation results later on...
    cnvData$Indiv <- paste0(cnvData$Indiv, ".hq.sorted.marked")
    cnvData$sample <- paste0(cnvData$sample, ".hq.sorted.marked")
    outputFile <- file.path(outputFolder, "all_cnv_calls.txt")
    write.table(cnvData, file = outputFile, sep='\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

    # GenomicRanges object----
    saveResultsFileToGR(outputFolder, "all_cnv_calls.txt", chrColumn = "CHROM", sampleColumn = "sample",
                       startColumn = "start", endColumn = "end", cnvTypeColumn = "copy_number")

    #Temporary files
    #Delete temporary files if specified
    if(keepTempFiles == "false"){
      filesAll <- list.files(outputFolder, full.names = TRUE)
      filesToKeep <- c("failedROIs.csv", "grPositives.rds", "cnvs_summary.tsv", "cnvFounds.csv", "cnvFounds.txt", "all_cnv_calls.txt", "calls_all.txt", "failures_Failures.txt", "cnv_calls.tsv")
      filesToRemove <- list(filesAll[!(filesAll %in% grep(paste(filesToKeep, collapse = "|"), filesAll, value = TRUE))])
      do.call(unlink, filesToRemove)
      foldersToRemove <- list.dirs(outputFolder, full.names = TRUE, recursive = FALSE)
      unlink(foldersToRemove, recursive = TRUE)
    }
  }
}


print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)
