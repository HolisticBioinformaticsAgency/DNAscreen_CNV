#!/usr/bin/env bash

## let's make sure it runs on one run of BAM folder before implementing a loop for multiple runs
## maybe multiple runs can be done on nextflow pipeline


# set -euo pipefail # use this when debugging...
set -x   # (prints every command before running it)

if [ $# -lt 1 ]; then
  echo "Usage: $0 clearCNV_params.yaml"
  exit 1
fi

CONFIG=$1

# --- Load YAML values ---
OUTPUT=$(yq -r '.outputFolder' $CONFIG)
PROJECT=$(yq -r '.project_name' $CONFIG)
BAMS_DIR=$(yq -r '.bams_dir' $CONFIG)
BED=$(yq -r '.bed_file' $CONFIG)
FASTA=$(yq -r '.fasta_file' $CONFIG)
BLACKLIST=$(yq -r '.blacklist' $CONFIG)

CORES=$(yq -r '.cores' $CONFIG)
EXPECTED=$(yq -r '.expected_artefacts' $CONFIG)
SCORE=$(yq -r '.sample_score_factor' $CONFIG)
ZSCALE=$(yq -r '.zscale' $CONFIG)
SIZE=$(yq -r '.size' $CONFIG)
DEL_CUTOFF=$(yq -r '.del_cutoff' $CONFIG)
DUP_CUTOFF=$(yq -r '.dup_cutoff' $CONFIG)
TRANS=$(yq -r '.trans_prob' $CONFIG)
GROUPS=$(yq -r '.minimum_group_sizes' $CONFIG)

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
mkdir -p $OUTPUT/bams
BAM_TXT=$OUTPUT/bams/all_bams.txt
ls $BAMS_DIR/*.bam > $BAM_TXT

# --- Run clearCNV ---
source activate mamba-env

clearCNV workflow_cnv_calling \
  -w $OUTPUT \
  -p $PROJECT \
  -r $FASTA \
  -b $BAM_TXT \
  -d $BED \
  -k $BLACKLIST \
  -c $CORES \
  --expected_artefacts $EXPECTED \
  --sample_score_factor $SCORE \
  --minimum_group_sizes 20 \
  --zscale $ZSCALE \
  --size $SIZE \
  --del_cutoff $DEL_CUTOFF \
  --dup_cutoff $DUP_CUTOFF \
  --trans_prob $TRANS

# --- Post-process results ---
RESULTS=$OUTPUT/$PROJECT/results/cnv_calls.tsv
FINAL=$OUTPUT/${PROJECT}_cnv_calls.tsv


echo "✅ clearCNV finished. Results: $FINAL"
