#!/bin/bash

## Script for calculating FRIP score 
## 
## Example usage:
## bamDir=. peakDir=. outDir=. sbatch --array=0-0 calculate_frip_score.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=64G

# Job name and output
#SBATCH -J frip_score
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=8

# Load modules
module load bedtools
module load samtools

# Define query files
queries=($(ls ${bamDir}/*.sorted.bam | xargs -n 1 basename | sed 's/.sorted.bam//g'))

# Run the thing
pwd; hostname; date

echo "Bedtools version: "$(bedtools --version)
echo "Samtools version: "$(samtools --version)

echo "Processing file: "${bamDir}/${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Calculating total reads and peak reads...")

# Set a variable equal to the total number of reads (from bam file)
totalReads=`samtools view -c ${bamDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam`

echo "Total number of reads: "$(echo $totalReads)

# Set a variable equal to the number of reads in peaks (from macs2 narrowpeak file)
peakReads=`bedtools sort -i ${peakDir}/${queries[$SLURM_ARRAY_TASK_ID]}_mergedpeaks.bed \
| bedtools merge -i - \
| bedtools intersect -u -nonamecheck -a ${bamDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam -b - -ubam \
| samtools view -c`

echo "Number of reads in peaks: "$(echo $peakReads)

# Calculate FRiP score
echo $(date +"[%b %d %H:%M:%S] Calculating FRIP...")

FRIP=`awk "BEGIN {print "${peakReads}"/"${totalReads}"}"`

echo "${queries[$SLURM_ARRAY_TASK_ID]}" > ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_FRIP.txt
echo "Total number of mapped reads: "$(echo $totalReads) >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_FRIP.txt
echo "Number of mapped reads in peaks: "$(echo $peakReads) >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_FRIP.txt
echo "Fraction of mapped reads in peaks (FRiP): "$(echo $FRIP) >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_FRIP.txt

echo $(date +"[%b %d %H:%M:%S] Done!")
