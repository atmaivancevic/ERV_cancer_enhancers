#!/bin/bash

## Script for computing bam read count table
##
## Example usage:
## outDir=. sbatch compute_bam_count_table.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=1-00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J bamCountTable
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# load modules
module load bedtools

# run the thing
pwd; hostname; date

echo "Starting bedtools..."
echo $(date +"[%b %d %H:%M:%S] Computing matrix of bam read counts...")

bedtools multicov \
-bams UT_K27ac_rep1.sorted.bam UT_K27ac_rep2.sorted.bam Cobi_K27ac_rep1.sorted.bam Cobi_K27ac_rep2.sorted.bam TNF_K27ac_rep1.sorted.bam TNF_K27ac_rep2.sorted.bam \
-bed ALL_K27ac_CnR_top20k_peaks.bed \
> $outDir/H3K27ac_CnR_bam_counts_within_peaks.tab

echo $(date +"[%b %d %H:%M:%S] Done!")
