#!/bin/bash

## Script for merging narrowpeak files from the two modes 
## 
## Example usage:
## inDir1=3_macs2_output_PE/tf inDir2=4_macs2_output_SE/tf outDir=5_macs2_output_merged/tf sbatch --array=0-0 merge_peak_files.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem=16G

# Job name and output
#SBATCH -J merge_peaks
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Load modules
module load bedtools

# Define query files
queries=($(ls ${inDir1}/*_peaks.narrowPeak | xargs -n 1 basename))

# Run the thing
pwd; hostname; date

echo "Bedtools version: "$(bedtools --version)

echo "Processing file: "${inDir1}/${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Merging peak files...")

cat ${inDir1}/${queries[$SLURM_ARRAY_TASK_ID]} ${inDir2}/${queries[$SLURM_ARRAY_TASK_ID]} \
| bedtools sort -i - \
| bedtools merge -i - -c 5 -o max 
| awk '{print $1 "\t" $2 "\t" $3 "\t" "'${queries[$SLURM_ARRAY_TASK_ID]%_peaks.narrowPeak}_peak'"NR "\t" $4 "\t" "."}' \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_peaks.narrowPeak}_mergedpeaks.bed

echo $(date +"[%b %d %H:%M:%S] Done!")
