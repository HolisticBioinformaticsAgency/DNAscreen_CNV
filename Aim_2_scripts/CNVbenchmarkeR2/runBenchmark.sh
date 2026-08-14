#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --nodes=2
#SBATCH --job-name=benchmark
#SBATCH --output=benchmark.out
#SBATCH --error=benchmark.err

module load R
module load perl
module load java/openjdk-21.0.1

Rscript runBenchmark.R