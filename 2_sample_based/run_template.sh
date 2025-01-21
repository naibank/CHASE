#!/bin/bash

#SBATCH --job-name MGRB
#SBATCH -N 1 -c 8
#SBATCH --mem=48G
#SBATCH --time=48:00:00

module load R/4.4.0

Rscript 1_getCHEvent.R
