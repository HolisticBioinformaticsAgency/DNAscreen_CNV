#!/bin/bash
#SBATCH --job-name=cnv_framework
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --output=slurm-%j.out

# Usage:
#   sbatch run_framework.slurm                                  # runs 4-35
#   sbatch run_framework.slurm --samplesheet my_runs.csv        # any extra Nextflow args
#   sbatch run_framework.slurm -resume                          # resume a previous run

source /usr/share/Modules/init/bash

# Don't `module load nextflow/21.10.6`: it also loads java/openjdk-1.14.02, whose
# JAVA_HOME points to a JVM that no longer exists (CAPSULE EXCEPTION). Use the
# module's launcher directly with the system Java, and the cached 22.10.4
# framework (the version used for previous runs).
SYSTEM_JAVA=$(readlink -f /usr/bin/java)
export JAVA_HOME=$(dirname "$(dirname "$SYSTEM_JAVA")")
export JAVA_CMD="$SYSTEM_JAVA"
export PATH=/usr/local/nextflow/21.10.6/bin:$PATH
export NXF_VER=22.10.4

# Submit from the pipeline directory (sbatch copies this script to a spool dir)
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

nextflow run main.nf "$@"
