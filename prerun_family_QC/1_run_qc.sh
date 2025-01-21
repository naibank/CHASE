#!/bin/bash

#SBATCH --job-name CHASE_pre
#SBATCH -N 1 -c 8
#SBATCH --mem=48G
#SBATCH --time=48:00:00

module load R/3.5.1
Rscript sampleCNV_QCs.R
