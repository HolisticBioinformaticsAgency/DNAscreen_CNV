#!/bin/bash
#SBATCH --job-name=target_hifi
#SBATCH --output=target_hifi_%j.out
#SBATCH --error=target_hifi_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

# --- Environment setup --- #
source ~/.bashrc
conda activate pacbio_target

umask 002

BATCH=$1
BIOSAMPLES=$2
CCSREADS=$3
TARGETS=$4
PROBES=${5:-None}

if [ -z $3 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
 echo -e "\nUsage: sbatch $(basename $0) <batch_name> <biosample_csv> <hifi_reads> <target_bed> [<probe_bed>]\n"
 exit 0
fi

mkdir -p "batches/${BATCH}/"

# --- Use mkdir as a lock to prevent multiple simultaneous jobs --- #
LOCKDIR="batches/${BATCH}/process_batch.lock"

if ! mkdir "$LOCKDIR" 2>/dev/null; then
    echo "Batch $BATCH is already being processed. Exiting."
    exit 1
fi

# Ensure the lock is removed on exit or error
trap "rm -rf '$LOCKDIR'" EXIT SIGINT SIGTERM ERR

# --- execute snakemake --- #
snakemake --reason \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    --configfile /home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/demultiplex_pb/workflow_HiFiTargetEnrichment/config.yaml \
    --config batch="${BATCH}" \
             biosamples="${BIOSAMPLES}" \
             ccsReads="${CCSREADS}" \
             targets="${TARGETS}" \
             probes="${PROBES}" \
             scripts=/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/demultiplex_pb/workflow_HiFiTargetEnrichment/scripts \
    --nolock \
    --local-cores 4 \
    --jobs 750 \
    --max-jobs-per-second 1 \
    --use-conda --conda-frontend conda \
    --use-singularity --singularity-args '--nv ' \
    --latency-wait 300 \
    --cluster-config /home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/demultiplex_pb/workflow_HiFiTargetEnrichment/cluster.yaml \
    --cluster "sbatch --partition={cluster.partition} \
                      --cpus-per-task={cluster.cpus} \
                      --output={cluster.out} {cluster.extra} " \
    --snakefile /home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/demultiplex_pb/workflow_HiFiTargetEnrichment/Snakefile
