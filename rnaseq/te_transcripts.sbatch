#!/bin/bash

## Script for running TEtranscripts
##
## Example usage:
## sbatch te_transcripts.sbatch

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --time=48:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J TEtranscripts
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load TEtranscripts
module load python/2.7.14/tetranscripts/2.1.4

# Run TEtranscripts for each comparison, eg:

# (1) Cobimetinib vs Untreated
TEtranscripts --format BAM \
--sortByPos --mode multi \
-i 100 --stranded reverse \
--DESeq -n DESeq_default \
--GTF gencode.v34.annotation.gtf \
--TE GRCh38_GENCODE_rmsk_TE.gtf \
--project Cobi_24hr_bothReps \
-t Cobi_24_1_S54.k100.bam Cobi_24_2_S55.k100.bam \
-c UT_24_1_S52.k100.bam UT_24_2_S53.k100.bam

# (2) TNF-alpha vs Untreated
TEtranscripts --format BAM \
--sortByPos --mode multi \
-i 100 --stranded reverse \
--DESeq -n DESeq_default \
--GTF gencode.v34.annotation.gtf \
--TE GRCh38_GENCODE_rmsk_TE.gtf \
--project TNF_24hr_bothReps \
-t TNF_24_1_S58.k100.bam TNF_24_2_S59.k100.bam \
-c UT_24_1_S52.k100.bam UT_24_2_S53.k100.bam

echo $(date +"[%b %d %H:%M:%S] Done")
