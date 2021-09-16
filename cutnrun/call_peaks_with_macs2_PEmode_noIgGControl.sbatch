#!/bin/bash

## Script for running macs2 on paired-end data (without input igg)
##
## Example usage:
## inDir=. outDir=. sbatch --array=0-0 call_peaks_with_macs2_PEmode_noIgGControl.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=12:00:00
#SBATCH --mem=64G

# Job name and output
#SBATCH -J macs2_call_peaks_noigg
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Load modules
module load python/2.7.14/MACS/2.1.1

# Define query files
queries=($(ls ${inDir}/*.bam | xargs -n 1 basename))

# Run the thing
pwd; hostname; date

echo "macs2 version: "$(macs2 --version)
echo "Target file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Starting macs2...")

macs2 callpeak \
--format BAMPE \
--treatment ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
--pvalue 0.01 \
--name ${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam} \
--outdir ${outDir} \
-g hs \
--SPMR -B \
--call-summits

echo $(date +"[%b %d %H:%M:%S] Done!")
