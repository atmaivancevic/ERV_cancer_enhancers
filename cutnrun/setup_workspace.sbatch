#!/bin/bash

## Script to set up directory layout 
##
## Example usage:
## projectDir=. sbatch setup_workspace.q

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=0:10:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J setup_workspace
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Making project directories...")

mkdir -p ${projectDir}/1_trimmed_fastq
mkdir -p ${projectDir}/1_trimmed_fastq/control
mkdir -p ${projectDir}/1_trimmed_fastq/histone
mkdir -p ${projectDir}/1_trimmed_fastq/tf

mkdir -p ${projectDir}/2_bams
mkdir -p ${projectDir}/2_bams/control
mkdir -p ${projectDir}/2_bams/histone
mkdir -p ${projectDir}/2_bams/tf

mkdir -p ${projectDir}/3_macs2_output_PE
mkdir -p ${projectDir}/3_macs2_output_PE/histone
mkdir -p ${projectDir}/3_macs2_output_PE/tf

mkdir -p ${projectDir}/4_macs2_output_SE
mkdir -p ${projectDir}/4_macs2_output_SE/histone
mkdir -p ${projectDir}/4_macs2_output_SE/tf

mkdir -p ${projectDir}/5_macs2_output_merged
mkdir -p ${projectDir}/5_macs2_output_merged/histone
mkdir -p ${projectDir}/5_macs2_output_merged/tf

mkdir -p ${projectDir}/6_bigwigs

#mkdir -p ${projectDir}/7_seacr_output
#mkdir -p ${projectDir}/7_seacr_output/control
#mkdir -p ${projectDir}/7_seacr_output/histone
#mkdir -p ${projectDir}/7_seacr_output/tf

mkdir -p ${projectDir}/8_signal_heatmaps

mkdir -p ${projectDir}/9_giggle_output

mkdir -p ${projectDir}/reports
mkdir -p ${projectDir}/reports/bam_fragment_size
mkdir -p ${projectDir}/reports/frip_scores
mkdir -p ${projectDir}/reports/raw_fastqc
mkdir -p ${projectDir}/reports/raw_multiqc
mkdir -p ${projectDir}/reports/trimmed_fastqc
mkdir -p ${projectDir}/reports/trimmed_multiqc

echo $(date +"[%b %d %H:%M:%S] Done!")
