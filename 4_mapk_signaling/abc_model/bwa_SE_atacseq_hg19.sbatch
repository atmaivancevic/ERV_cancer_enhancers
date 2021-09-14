#!/bin/bash

## Script for aligning single-ended atacseq fastq files to hg19
##
## Example usage:
## inDir=. outDir=. sbatch --array 0-0 bwa_SE_atacseq_hg19.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=6:00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J bwa-single
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# define key variables
bwaIndexDir=/Shares/CL_Shared/db/genomes/hg19/index/bwa
bwaIndex=hg19.main.fa

# define query files
queries=($(ls $inDir/*trimmed.fastq.gz | xargs -n 1 basename))

# load modules
module load bwa samtools

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Starting bwa alignment...")

bwa mem -t 8 \
${bwaIndexDir}/${bwaIndex} \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
| samtools view -Sb -q 10 - \
| samtools sort -@ 8 - -o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}.sorted.uniq.bam

echo $(date +"[%b %d %H:%M:%S] Done!")
