#!/bin/bash

#SBATCH -N 1 -c 2
#SBATCH --mem=150G
#SBATCH --tmp=150G
#SBATCH -t 180:00:00

module load R/4.4.0

Rscript 2_get_CNV_SNV_table.R
