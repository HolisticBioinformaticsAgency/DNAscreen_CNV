# Generates summary file
#USAGE: Rscript summary.r [tools_params_file] [datasets_params_file]
source(if (basename(getwd()) == "utils") "cnvStats.r" else "utils/cnvStats.r") # Load class definitions
suppressPackageStartupMessages(library(yaml))

# load tools and datasets from yaml files
args <- commandArgs(TRUE)
cat("\n"); print(args)
if(length(args)>0) {
  toolsParamsFile <- args[1]
  datasetsParamsFile <- args[2]
} else {
  toolsParamsFile <- "tools.yaml"
  datasetsParamsFile <- "datasets.yaml"
}


#Load the parameters file
tools <- yaml.load_file(toolsParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)

# Run summary for all datasets and tools
for (dName in names(datasets)) {
  dataset <- datasets[[dName]]

  if (dataset$include){

    # Create SummaryStats object and load current dataset validated results
    ss <- SummaryStats()
    ss$loadValidatedResults(path = dataset$validated_results_file,
                            datasetName = dName,
                            bedFile = dataset$bed_file_summary)

    # Calculate metrics
    for (alg in names(tools)) {
      useAlg <- tools[[alg]]
      if (useAlg){
        outputFolder <- file.path(getwd(), "output", paste0(alg, "-", dName))
        ss$loadAlgorithmResults(outputFolder = outputFolder,
                                algorithmName = alg,
                                datasetName = dName,
                                bedFile = dataset$bed_file_summary)
      }
    }

    # Write results to summary file
    ss$writeSummary(file.path(getwd(), "output", "summary", paste0("summary-", dName, ".txt")), dName)
    ss$writeCSVresults(file.path(getwd(), "output", "summary", paste0("results-", dName, ".tsv")), dName)
  }
}
