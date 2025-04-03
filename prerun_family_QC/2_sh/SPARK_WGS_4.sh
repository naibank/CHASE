#!/bin/bash

#SBATCH --job-name SPARK_WGS_4
#SBATCH -N 1 -c 8
#SBATCH --mem=48G
#SBATCH --time=48:00:00

module load R/3.5.1

Rscript filterMosaicSVs.R SPARK_WGS_4
