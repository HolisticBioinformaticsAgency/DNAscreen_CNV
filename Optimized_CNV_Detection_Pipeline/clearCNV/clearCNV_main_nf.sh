#!/usr/bin/env bash

# Exit on error, undefined vars, and pipe failures
# set -euo pipefail
set -x   # print each command before running it

if [ $# -lt 1 ]; then
  echo "Usage: $0 clearCNV_params_nf.yaml [--bams_dir DIR] [--out DIR]"
  exit 1
fi

CONFIG=$1
shift

# --- Load YAML values ---
PROJECT=$(yq -r '.project_name' "$CONFIG")
BED=$(yq -r '.bed_file' "$CONFIG")
FASTA=$(yq -r '.fasta_file' "$CONFIG")
BLACKLIST=$(yq -r '.blacklist' "$CONFIG")
CORES=$(yq -r '.cores' "$CONFIG")
EXPECTED=$(yq -r '.expected_artefacts' "$CONFIG")
SCORE=$(yq -r '.sample_score_factor' "$CONFIG")
ZSCALE=$(yq -r '.zscale' "$CONFIG")
SIZE=$(yq -r '.size' "$CONFIG")
DEL_CUTOFF=$(yq -r '.del_cutoff' "$CONFIG")
DUP_CUTOFF=$(yq -r '.dup_cutoff' "$CONFIG")
TRANS=$(yq -r '.trans_prob' "$CONFIG")
GROUPS=$(yq -r '.minimum_group_sizes' "$CONFIG")

# --- Allow overrides from CLI ---
BAMS_DIR=""
OUTPUT=""
while [[ $# -gt 0 ]]; do
  case $1 in
    --bams_dir) BAMS_DIR="$2"; shift 2;;
    --out) OUTPUT="$2"; shift 2;;
    *) echo "Unknown option: $1"; exit 1;;
  esac
done

# --- Set defaults if not provided ---
BAMS_DIR=${BAMS_DIR:-"./bams"}
OUTPUT=${OUTPUT:-"./output"}

echo "Parsed parameters:"
echo "OUTPUT=$OUTPUT"
echo "PROJECT=$PROJECT"
echo "BAMS_DIR=$BAMS_DIR"
echo "BED=$BED"
echo "FASTA=$FASTA"
echo "BLACKLIST=$BLACKLIST"
echo "CORES=$CORES"
echo "EXPECTED=$EXPECTED"
echo "SCORE=$SCORE"
echo "GROUPS=$GROUPS"
echo "ZSCALE=$ZSCALE"
echo "SIZE=$SIZE"
echo "DEL_CUTOFF=$DEL_CUTOFF"
echo "DUP_CUTOFF=$DUP_CUTOFF"
echo "TRANS=$TRANS"

# --- Prepare output dirs ---
mkdir -p "$OUTPUT/bams"
BAM_TXT="$OUTPUT/bams/all_bams.txt"

# --- Check for BAM files ---
shopt -s nullglob
bam_files=("$BAMS_DIR"/*.bam)
if [ ${#bam_files[@]} -eq 0 ]; then
    echo "No BAM files found in $BAMS_DIR"
    exit 1
fi

printf "%s\n" "${bam_files[@]}" > "$BAM_TXT"

# --- Activate environment ---
source activate mamba-env

# --- Run clearCNV ---
clearCNV workflow_cnv_calling \
  -w "$OUTPUT" \
  -p "$PROJECT" \
  -r "$FASTA" \
  -b "$BAM_TXT" \
  -d "$BED" \
  -k "$BLACKLIST" \
  -c "$CORES" \
  --expected_artefacts "$EXPECTED" \
  --sample_score_factor "$SCORE" \
  --minimum_group_sizes "$GROUPS" \
  --zscale "$ZSCALE" \
  --size "$SIZE" \
  --del_cutoff "$DEL_CUTOFF" \
  --dup_cutoff "$DUP_CUTOFF" \
  --trans_prob "$TRANS"

# --- Post-process results ---
RESULTS="$OUTPUT/$PROJECT/results/cnv_calls.tsv"
FINAL="$OUTPUT/cnv_calls.tsv"

awk 'BEGIN{FS=OFS="\t"} NR==1{print $0,"CNV.type"; next}
     {if($5=="DEL") t="deletion"; else t="duplication"; print $0,t}' \
     "$RESULTS" > "$FINAL"

echo "✅ clearCNV finished. Results: $FINAL"
