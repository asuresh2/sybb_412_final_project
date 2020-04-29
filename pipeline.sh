#!/bin/bash
#SBATCH -J pipeline                  # A single job name for the array
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=64000                 # Memory request (4Gb)
#SBATCH -t 1-00:00                  # Maximum execution time (D-HH:MM)


module load gcc/7.3.0
module load java
module load python

# Set up path variables.
PATH=$PATH:/mnt/pan/courses/sybb412/tools/FastQC/
PATH=$PATH:/home/axs1114/.local/bin/
PATH=$PATH:/mnt/pan/courses/sybb412/axs1114/TrimGalore-0.6.5
PATH=$PATH:/mnt/pan/courses/sybb412/axs1114/final_project/scripts/subread-2.0.0-Linux-x86_64/bin/
PATH=$PATH:/mnt/pan/courses/sybb412/axs1114/final_project/scripts/STAR-2.7.3a/bin/Linux_x86_64/

# Specify data directory.
data_dir="/mnt/pan/courses/sybb412/axs1114/final_project/data"

# Specify sample name.
sample_name=$1
echo ${sample_name}


# Make output directories for fastqc output and alignment output.
mkdir /mnt/pan/courses/sybb412/axs1114/final_project/data/fastqc_output/${sample_name}
fastqc_output_dir="/mnt/pan/courses/sybb412/axs1114/final_project/data/fastqc_output/${sample_name}"
bam_file_output="/mnt/pan/courses/sybb412/axs1114/final_project/data/bam_files"
sorted_bam_files="/mnt/pan/courses/sybb412/axs1114/final_project/data/bam_files/bam_sorted"

# Run fastqc and trimming.
trim_galore -j 4 --paired --fastqc ${data_dir}/fastq_files/${sample_name}/*_1.fastq ${data_dir}/fastq_files/${sample_name}/*_2.fastq --output_dir ${fastqc_output_dir}  #--path_to_cutadapt /home/axs1114/.local/bin/cutadapt/



# Find index file for alignment.

module load gcc/6.3.0
module load samtools
module load hisat2

hisat2 -x /mnt/pan/courses/sybb412/axs1114/final_project/data/reference_genomes/grch38_tran/genome_tran -p 4 -1 ${fastqc_output_dir}/${sample_name}_1_val_1.fq -2 ${fastqc_output_dir}/${sample_name}_2_val_2.fq | samtools view -bS | samtools sort -n - > /mnt/pan/courses/sybb412/axs1114/final_project/data/bam_files/${sample_name}.bam


gtf_file=/mnt/pan/courses/sybb412/axs1114/final_project/data/hg38.ncbiRefSeq.gtf

samtools sort -n ${bam_file_output}/${sample_name}.bam -o ${bam_file_output}/bam_sorted/${sample_name}_sorted.bam

featureCounts -p -t exon -g gene_id -M -O -a $gtf_file -T 4 -s 1 -o /mnt/pan/courses/sybb412/axs1114/final_project/data/counts_sorted_HISAT2_refseq.txt ${sorted_bam_files}/*.bam


## If STAR workflow is desired, use the following commands for quantification of transcripts. ##

#star_output_files="/mnt/pan/courses/sybb412/axs1114/final_project/data/star_output"

#mkdir /mnt/pan/courses/sybb412/axs1114/final_project/data/STAR_index

# Generate STAR index
#STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /mnt/pan/courses/sybb412/axs1114/final_project/data/STAR_index --genomeFastaFiles /mnt/pan/courses/sybb412/axs1114/final_project/data/reference_genomes/GRCh38.p13.genome.fa --sjdbGTFfile /mnt/pan/courses/sybb412/axs1114/final_project/data/chess2.2.gtf

# Run alignment with STAR.
#STAR --runThreadN 4 --genomeDir /mnt/pan/courses/sybb412/axs1114/final_project/data/STAR_index --readFilesIn ${fastqc_output_dir}/${sample_name}_1_val_1.fq ${fastqc_output_dir}/${sample_name}_2_val_2.fq --outFileNamePrefix ${star_output_files}/${sample_name} --outSAMtype BAM Unsorted




