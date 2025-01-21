#!/bin/bash

#SBATCH -N 1 -c 2
#SBATCH --mem=100G
#SBATCH --tmp=100G
#SBATCH -t 180:00:00

root_dir=/hpf/largeprojects/tcagstor/scratch/kara.han/CHASE
out_dir=$root_dir/data/SNV
tool=$root_dir/script/1_get_SNVs.R

# /hpf/largeprojects/tcagstor/tools/data/MSSNG/CG/variants/SNVs+indels/exonic+splicing 1738
# /hpf/largeprojects/tcagstor/tools/data/MSSNG/ILMN/variants/SNVs+indels/exonic+splicing 9618
# /hpf/largeprojects/tcagstor/tools/data/SSC/variants/SNVs+indels/exonic+splicing 9205
# /hpf/largeprojects/tcagstor/tools/data/SPARK_WGS_1/variants/SNVs+indels/exonic+splicing 2629
# /hpf/largeprojects/tcagstor/tools/data/SPARK_WGS_2/variants/SNVs+indels/exonic+splicing 2365
# /hpf/largeprojects/tcagstor/tools/data/SPARK_WGS_3/variants/SNVs+indels/exonic+splicing 2867
# /hpf/largeprojects/tcagstor/tools/data/SPARK_WGS_4/variants/SNVs+indels/exonic+splicing 3684
infile_dir=/hpf/largeprojects/tcagstor/tools/data/MSSNG/CG/variants/SNVs+indels/exonic+splicing
cohort=MSSNG_CG
cohort_outdir=$out_dir/$cohort
mkdir $out_dir
mkdir $cohort_outdir

# Divide the total num of files into chunks
num_smpl=$(find $infile_dir -name '*.tsv.gz' | wc -l)
echo "There are $num_smpl samples"
num_per_chunk=150
a=$(($num_smpl/$num_per_chunk))
b=$(($num_smpl%$num_per_chunk))
echo "a is $a"
echo "b is $b"
i=1
while [ "$i" -le "$a" ];
do
    echo "i is $i"
    c=$(($i*$num_per_chunk))
    files=$(find $infile_dir -name '*.tsv.gz' | head -n $c | tail -n $num_per_chunk) 
    snv_subset_outpath=$cohort_outdir/"$cohort"_snv_subset$i.txt
    sbatch --export=tool=$tool,files="${files}",snv_subset_outpath=$snv_subset_outpath $root_dir/script/1_get_SNVs.sh
    i=$(($i+1))
done
echo "i is $i"
files=$(find $infile_dir -name '*.tsv.gz' | tail -n $b)
snv_subset_outpath=$cohort_outdir/"$cohort"_snv_subset$i.txt
sbatch --export=tool=$tool,files="${files}",snv_subset_outpath=$snv_subset_outpath $root_dir/script/1_get_SNVs.sh
echo "Finished"
