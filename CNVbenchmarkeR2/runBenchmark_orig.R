# Description: Run the benchmark analysis
# USAGE: Rscript runBenchmark.R [-t tools_file] [-d datasets_file] [-f keepTempFiles]
# keepTempFiles: if true, temp files will not be removed (Default: true)
# libs
library("optparse")
library(yaml)
start.benchmark <- format(Sys.time(), "%X_%Y")
# Check runBenchmark.R is being called from CNVbenchmarkeR2 folder
if (length(list.files(pattern = "runBenchmark.R"))== 0){
  cat("Sorry, runBenchmark.R should be called from CNVbenchmarkeR2 folder\n")
  quit()
}

# Build options list
option_list <- list(
  make_option(c("-t", "--tools"), type="character", default="tools.yaml",
              help="Path to tools file (yaml)", metavar="character"),
  make_option(c("-d", "--datasets"), type="character", default="datasets.yaml",
              help="Path to datasets file (yaml)", metavar="character"),
  make_option(c("-f", "--keepTempFiles"), type="character", default="true",
              help="Include temporary files in the output folders", metavar="character")
);
opt_parser <- OptionParser(option_list=option_list);

# Load params
args <- parse_args(opt_parser);
tools <- yaml.load_file(args$tools)



#create logs/output folder if not exists
dir.create("logs", showWarnings = F)
dir.create("output", showWarnings = F)


## Execute tools on selected datasets ##

# AtlasCNV
if (tools$atlasCNV == T){
  cat(as.character(Sys.time()), " - Executing Atlas-CNV\n")
  cmd <- paste("Rscript tools/atlasCNV/runAtlascnv.r tools/atlasCNV/atlasCNVParams.yaml", args$datasets, args$keepTempFiles, "false", ">",  paste0("logs/atlasCNV_", start.benchmark, ".log 2>&1"))
  system(cmd)
}

# ClearCNV
if (tools$clearCNV == T){
  cat(as.character(Sys.time()), " - Executing clearCNV\n")
  cmd <- paste("Rscript tools/clearCNV/runClearCNV.r tools/clearCNV/clearCNVParams.yaml", args$datasets, args$keepTempFiles," >",  paste0("logs/clearCNV_", start.benchmark, ".log 2>&1"))
  system(cmd)
}

# ClinCNV
if (tools$clincnv == T){
  cat(as.character(Sys.time()), " - Executing ClinCNV\n")
  cmd <- paste("Rscript tools/clincnv/runClincnv.R tools/clincnv/clincnvParams.yaml", args$datasets, args$keepTempFiles, " >",  paste0("logs/clincnv_", start.benchmark, ".log 2>&1"))
  system(cmd)
}

# Cnvkit
if (tools$cnvkit == T){
  cat(as.character(Sys.time()), " - Executing cnvkit\n")
  cmd <- paste("Rscript tools/cnvkit/runCnvkit.r tools/cnvkit/cnvkitParams.yaml", args$datasets, args$keepTempFiles," >",  paste0("logs/cnvkit_", start.benchmark, ".log 2>&1"))
  system(cmd)
}

# Cobalt
if (tools$cobalt == T){
  cat(as.character(Sys.time()), " - Executing Cobalt\n")
  cmd <- paste("Rscript tools/cobalt/runCobalt.r tools/cobalt/cobaltParams.yaml", args$datasets, args$keepTempFiles," >",  paste0("logs/cobalt_", start.benchmark, ".log 2>&1"))
  system(cmd)
}

# CODEX2
if (tools$codex2 == T){
  cat(as.character(Sys.time()), " - Executing CODEX2\n")
  cmd <- paste("Rscript tools/codex2/runCodex2.r tools/codex2/codex2Params.yaml", args$datasets, args$keepTempFiles," >",  paste0("logs/codex2_", start.benchmark, ".log 2>&1"))
  system(cmd)
}

# CoNVaDING
if (tools$convading == T){
  cat(as.character(Sys.time()), " - Executing CoNVaDING\n")
  cmd <- paste("Rscript tools/convading/runConvading.r tools/convading/convadingParams.yaml", args$datasets, args$keepTempFiles," >",  paste0("logs/convading_", start.benchmark, ".log 2>&1"))
  system(cmd)
}

# Decon
if (tools$decon == T){
  cat(as.character(Sys.time()), " - Executing DECoN\n")
  cmd <- paste("Rscript tools/decon/runDecon.r tools/decon/deconParams.yaml", args$datasets, args$keepTempFiles, " >",  paste0("logs/decon_", start.benchmark, ".log 2>&1"))
  system(cmd)
}

# ExomeDepth
if (tools$exomedepth == T){
  cat(as.character(Sys.time()), " - Executing ExomeDepth\n")
  cmd <- paste("Rscript tools/exomedepth/runExomedepth.r tools/exomedepth/exomedepthParams.yaml", args$datasets, args$keepTempFiles," >",  paste0("logs/exomedepth_", start.benchmark, ".log 2>&1"))
  system(cmd)
}

# GATK-gCNVcaller
if (tools$germlineCNVcaller == T){
  cat(as.character(Sys.time()), " - Executing germlineCNVcaller\n")
  cmd <- paste("Rscript tools/germlineCNVcaller/runGermlineCNVcaller.r tools/germlineCNVcaller/germlineCNVcallerParams.yaml", args$datasets, args$keepTempFiles,"false"," >",  paste0("logs/germlineCNVcaller_", start.benchmark, ".log 2>&1"))
  system(cmd)
}

# Panelcn.mops
if (tools$panelcnmops == T){
  cat(as.character(Sys.time()), " - Executing panelcn.MOPS\n")
  cmd <- paste("Rscript tools/panelcnmops/runPanelcnmops.r tools/panelcnmops/panelcnmopsParams.yaml", args$datasets, args$keepTempFiles, " >",  paste0("logs/panelcnmops_", start.benchmark, ".log 2>&1"))
  system(cmd)
}

# Viscap
if (tools$viscap == T){
  cat(as.character(Sys.time()), " - Executing Viscap\n")
  cmd <- paste("Rscript tools/viscap/runVisCap.r tools/viscap/viscapParams.yaml", args$datasets, args$keepTempFiles, "false", ">",  paste0("logs/viscap_", start.benchmark, ".log 2>&1"))
  system(cmd)
}

##  Generate summary file  ##

cat(as.character(Sys.time()), " - Generating summary file")
cmd <- paste("Rscript utils/summary.r", args$tools, args$datasets, ">",  paste0("logs/summary_", start.benchmark, ".log 2>&1"))
system(cmd)