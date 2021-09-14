#!/bin/bash

## Script for running bwa
##
## Example usage:
## inDir=. outDir=. sbatch --array 0-0 bwa_PE_cutnrun_hg19.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=8:00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J bwa_PE
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Define key variables
bwaIndexDir=/Shares/CL_Shared/db/genomes/hg19/index/bwa
bwaIndex=hg19.main.fa

# Set constant variables
numThreads=8

# Define query files
queries=($(ls ${inDir}/*fastq.gz | xargs -n 1 basename | sed 's/_R1_trimmed.fastq.gz//g' | sed 's/_R2_trimmed.fastq.gz//g' | uniq))

# Load modules
module load bwa
module load samtools
module load bedtools

# Print time and date
pwd; hostname; date

# Print program versions
echo "bwa version: "$(bwa)
echo "Samtools version: "$(samtools --version)
echo "Bedtools version: "$(bedtools --version)

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Starting bwa alignment...")

# Run bwa and filter with samtools, saving output bam
bwa mem \
-t ${numThreads} \
${bwaIndexDir}/${bwaIndex} \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R1_trimmed.fastq.gz \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R2_trimmed.fastq.gz \
| samtools view -Sb -q 10 -F 4 - \
| samtools sort -@ ${numThreads} - -o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}.sorted.bam

# Note that this will filter out low quality reads (-q 10), unmapped reads (-F 4), and chrM (samtools view -b - ${nonChrM}), then sort the bam.

# Create a bai index file
samtools index ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Done!")
