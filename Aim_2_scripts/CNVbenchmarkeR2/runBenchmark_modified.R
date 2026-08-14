# Description: Run the benchmark analysis
# USAGE: Rscript runBenchmark.R [tools_file] [datasets_file] [keepTempFiles]
# keepTempFiles: if true, temp files will not be removed (Default: true)

library(yaml)

start.benchmark <- format(Sys.time(), "%X_%Y")

# Ensure script is run from the correct directory
if (length(list.files(pattern = "runBenchmark.R")) == 0) {
  cat("Sorry, runBenchmark.R should be called from CNVbenchmarkeR2 folder\n")
  quit()
}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
tools_file <- ifelse(length(args) >= 1, args[1], "tools.yaml")
datasets_file <- ifelse(length(args) >= 2, args[2], "datasets.yaml")
keepTempFiles <- ifelse(length(args) >= 3, args[3], "true")

# Load tools configuration
tools <- yaml.load_file(tools_file)

# Create logs/output folders if not exist
dir.create("logs", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

# Function to run a tool
run_tool <- function(tool_name, script_path, params_file) {
  if (tools[[tool_name]] == TRUE) {
    cat(as.character(Sys.time()), " - Executing", tool_name, "\n")
    
    # Override keepTempFiles for specific tools
    tool_keepTempFiles <- ifelse(tool_name %in% c("germlineCNVcaller", "viscap"), "true", keepTempFiles)
    
    cmd <- paste("Rscript", script_path, params_file, , tool_keepTempFiles, ">", 
                 paste0("logs/", tool_name, "_", start.benchmark, ".log 2>&1"))
    system(cmd)
  }
}


# Execute tools
run_tool("atlasCNV", "tools/atlasCNV/runAtlascnv.r", "tools/atlasCNV/atlasCNVParams.yaml")
run_tool("clearCNV", "tools/clearCNV/runClearCNV.r", "tools/clearCNV/clearCNVParams.yaml")
run_tool("clincnv", "tools/clincnv/runClincnv.R", "tools/clincnv/clincnvParams.yaml")
run_tool("cnvkit", "tools/cnvkit/runCnvkit.r", "tools/cnvkit/cnvkitParams.yaml")
run_tool("cobalt", "tools/cobalt/runCobalt.r", "tools/cobalt/cobaltParams.yaml")
run_tool("codex2", "tools/codex2/runCodex2.r", "tools/codex2/codex2Params.yaml")
run_tool("convading", "tools/convading/runConvading.r", "tools/convading/convadingParams.yaml")
run_tool("decon", "tools/decon/runDecon.r", "tools/decon/deconParams.yaml")
run_tool("exomedepth", "tools/exomedepth/runExomedepth.r", "tools/exomedepth/exomedepthParams.yaml")
run_tool("germlineCNVcaller", "tools/germlineCNVcaller/runGermlineCNVcaller.r", "tools/germlineCNVcaller/germlineCNVcallerParams.yaml")
run_tool("panelcnmops", "tools/panelcnmops/runPanelcnmops.r", "tools/panelcnmops/panelcnmopsParams.yaml")
run_tool("viscap", "tools/viscap/runVisCap.r", "tools/viscap/viscapParams.yaml")

# Generate summary file
cat(as.character(Sys.time()), " - Generating summary file\n")
system(paste("Rscript utils/summary.r", tools_file, datasets_file, ">", paste0("logs/summary_", start.benchmark, ".log 2>&1")))
