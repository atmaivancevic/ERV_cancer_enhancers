#!/bin/bash

## Script for running bwa
##
## Example usage:
## inDir=. outDir=. sbatch --array 0-0 bwa_PE.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --time=8:00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J bwa_PE
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=4
nonChrM=$(cat ${genomeChrFile} | awk '{print $1}' | grep -v chrM | tr '\n' ' ')

# Define key variables
bwaIndexDir=/Shares/CL_Shared/db/genomes/hg38/index/bwa
bwaIndex=hg38.main.fa
genomeChrFile=/Shares/CL_Shared/db/genomes/hg38/fa/hg38.main.chrom.sizes

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
| samtools view -b - ${nonChrM} \
| samtools sort -@ ${numThreads} - \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}.sorted.bam

# Note that this will filter out low quality reads (-q 10), unmapped reads (-F 4), and chrM (samtools view -b - ${nonChrM}), then sort the bam.

# Create a bai index file
samtools index ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Done!")
