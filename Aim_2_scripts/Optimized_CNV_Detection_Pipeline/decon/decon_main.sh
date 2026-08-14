#!/usr/bin/env bash

module load R

# Runs DECoN using a config YAML file
# USAGE: ./run_decon.sh decon_params.yaml

set -x  # print commands for debugging
set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 decon_params.yaml"
  exit 1
fi

CONFIG=$1

# --- Load YAML values ---
OUTPUT=$(yq -r '.outputFolder' $CONFIG)
PROJECT=$(yq -r '.project_name' $CONFIG)
BAMS_DIR=$(yq -r '.bams_dir' $CONFIG)
BED=$(yq -r '.bed_file' $CONFIG)
FASTA=$(yq -r '.fasta_file' $CONFIG)
DECON_FOLDER=$(yq -r '.deconFolder' $CONFIG)
MINCORR=$(yq -r '.mincorr' $CONFIG)
MINCOV=$(yq -r '.mincov' $CONFIG)
TRANSPROB=$(yq -r '.transProb' $CONFIG)
KEEP_TEMP=$(yq -r '.keepTempFiles // "true"' $CONFIG)

echo "Parsed parameters:"
echo "OUTPUT=$OUTPUT"
echo "PROJECT=$PROJECT"
echo "BAMS_DIR=$BAMS_DIR"
echo "BED=$BED"
echo "FASTA=$FASTA"
echo "DECON_FOLDER=$DECON_FOLDER"
echo "MINCORR=$MINCORR"
echo "MINCOV=$MINCOV"
echo "TRANSPROB=$TRANSPROB"
echo "KEEP_TEMP=$KEEP_TEMP"

# --- Prepare output dirs ---
mkdir -p "$OUTPUT/$PROJECT/bams"
BAM_TXT="$OUTPUT/$PROJECT/bams/all_bams.txt"
find "$BAMS_DIR" -maxdepth 1 -name '*.bam' > "$BAM_TXT"

if [ ! -s "$BAM_TXT" ]; then
    echo "⚠️ No BAM files found in $BAMS_DIR"
    exit 1
fi

# --- DECoN workflow ---
echo "Starting DECoN for project $PROJECT at $(date)"

OUTPUT_BAMS="$OUTPUT/$PROJECT/output.bams"
OUTPUT_RDATA="$OUTPUT/$PROJECT/output.bams.RData"

# Pre-calc phase
echo "Running ReadInBams.R..."
Rscript "$DECON_FOLDER/ReadInBams.R" --bams "$BAMS_DIR" --bed "$BED" --fasta "$FASTA" --out "$OUTPUT_BAMS"
echo "ReadInBams.R finished"

# Identify failures
echo "Running IdentifyFailures.R..."
Rscript "$DECON_FOLDER/IdentifyFailures.R" --RData "$OUTPUT_RDATA" --mincorr "$MINCORR" --mincov "$MINCOV" --out "$OUTPUT/$PROJECT/failures"
echo "IdentifyFailures.R finished"

# Make CNV calls
echo "Running makeCNVcalls.R..."
Rscript "$DECON_FOLDER/makeCNVcalls.R" --RData "$OUTPUT_RDATA" --transProb "$TRANSPROB" --plot None --out "$OUTPUT/$PROJECT" --failures "$OUTPUT/$PROJECT/failures_Failures.txt"
echo "makeCNVcalls.R finished"

# Save CNV calls
## script works, but the calls_all.txt file is not in the calls directory.
RESULTS="$OUTPUT/${PROJECT}/calls_all.txt"
FINAL="$OUTPUT/${PROJECT}/${PROJECT}_cnv_calls.tsv"

if [ -f "$RESULTS" ]; then
    cp "$RESULTS" "$FINAL"
    echo "✅ DECoN finished. CNV calls saved to $FINAL"
else
    echo "⚠️ CNV calls not found!"
fi

# --- Cleanup temporary files ---
if [ "$KEEP_TEMP" = "false" ]; then
    echo "Removing temporary files..."
    find "$OUTPUT" -type f ! -name "${PROJECT}_cnv_calls.tsv" -delete
fi

echo "DECoN finished at $(date)"
