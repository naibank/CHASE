#!/bin/bash

#SBATCH -N 1 -c 2
#SBATCH --mem=100G
#SBATCH --tmp=100G
#SBATCH -t 180:00:00

root_dir="path_to_root_dir"
out_dir="path_to_output_directory"
tool="path_to_1_get_SNVs.R"

infile_dir="path_to_SNV_directory"
cohort="cohort_name"
cohort_outdir=$out_dir/$cohort
mkdir $out_dir
mkdir $cohort_outdir

# Divide the total num of files into chunks
num_smpl=$(find $infile_dir -name '*.tsv.gz' | wc -l) # Change *.tsv.gz to corresponding file suffix
echo "There are $num_smpl samples"
num_per_chunk=150 # Modify according to total number of samples
a=$(($num_smpl/$num_per_chunk))
b=$(($num_smpl%$num_per_chunk))
echo "a is $a"
echo "b is $b"
i=1
while [ "$i" -le "$a" ];
do
    echo "i is $i"
    c=$(($i*$num_per_chunk))
    files=$(find $infile_dir -name '*.tsv.gz' | head -n $c | tail -n $num_per_chunk) # Change *.tsv.gz to corresponding file suffix
    snv_subset_outpath=$cohort_outdir/"$cohort"_snv_subset$i.txt
    sbatch --export=tool=$tool,files="${files}",snv_subset_outpath=$snv_subset_outpath "path_to_1_get_SNVs.sh"
    i=$(($i+1))
done
echo "i is $i"
files=$(find $infile_dir -name '*.tsv.gz' | tail -n $b)
snv_subset_outpath=$cohort_outdir/"$cohort"_snv_subset$i.txt
sbatch --export=tool=$tool,files="${files}",snv_subset_outpath=$snv_subset_outpath "path_to_1_get_SNVs.sh"
echo "Finished"
