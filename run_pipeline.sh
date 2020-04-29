#!/bin/sh

sample_list="/mnt/pan/courses/sybb412/axs1114/final_project/scripts/SRR_Acc_List_last.txt"

cat $sample_list | while read i
do
	echo $i
	sbatch pipeline.sh ${i}
done