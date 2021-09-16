#!/bin/bash

## Script for running deeptools 
##
## Example usage:
## inDir=. outDir=. sbatch --array 0-0 get_fragment_size.q

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J get_fragment_size
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

pwd; hostname; date

# Define query files
queries=($(ls ${inDir}/*.sorted.bam | xargs -n 1 basename))

# Load modules
module load singularity

# Define key variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img
numThreads=8

# Get fragment size 
echo $(date +"[%b %d %H:%M:%S] Starting deeptools bamPEFragmentSize...")

singularity exec --bind /scratch/Users ${deeptools} \
bamPEFragmentSize \
--bamfiles ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
--table ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}_metrics.txt \
--histogram ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}_hist.png \
--samplesLabel ${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam} \
--plotTitle "Fragment size of PE sorted bam" \
--outRawFragmentLengths ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}_rawFragLengths.tab \
--numberOfProcessors ${numThreads} 

echo $(date +"[%b %d %H:%M:%S] Done!")
