#!/bin/bash

#SBATCH -N 1 -c 2
#SBATCH --mem=100G
#SBATCH --tmp=100G
#SBATCH -t 180:00:00

module load R/4.4.0

for file in $files
do
    echo "Processing file $file"
    Rscript $tool $file $snv_subset_outpath
done
