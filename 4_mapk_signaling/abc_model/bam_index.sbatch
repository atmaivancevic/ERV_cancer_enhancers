#!/bin/bash

## Script for aligning single-ended fastq files
##
## Example usage:
## sbatch --array 0-0 bam_index.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=6:00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J bam_index
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# define query files
queries=($(ls *.sorted.uniq.bam | xargs -n 1 basename))

# load modules
module load samtools

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Starting samtools index...")

# Create a bai index file
samtools index ${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done!")
