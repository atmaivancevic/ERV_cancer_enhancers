#!/bin/bash

## Script for subsetting bams by fragment size
##
## Example usage:
## inDir=. outDir=. sbatch --array 0-3 subset_by_fragment_size.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=1:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J subset_bam_by_size
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

pwd; hostname; date

# Define query files
queries=($(ls ${inDir}/*.sorted.bam | xargs -n 1 basename))

# Load modules
module load singularity
module load samtools

# define key variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img
numThreads=8

# Use "alignmentSieve" to filter reads by size (<120 or >120bp)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Starting deeptools alignmentSieve...")

singularity exec --bind /scratch/Users ${deeptools} \
alignmentSieve \
--bam ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
--maxFragmentLength 150 \
--outFile ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.lessThan150.sorted.bam.tmp \
--filterMetrics ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.lessThan150_metrics.txt \
--filteredOutReads ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.moreThan150.sorted.bam.tmp \
--numberOfProcessors ${numThreads}

# Sort and index the output bams
echo $(date +"[%b %d %H:%M:%S] Sorting bams (less than 150)...")

samtools sort \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.lessThan150.sorted.bam.tmp \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.lessThan150.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Removing temporary files...")

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.lessThan150.sorted.bam.tmp
rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.moreThan150.sorted.bam.tmp

echo $(date +"[%b %d %H:%M:%S] Indexing bams (less than 150)...")

samtools index \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.lessThan150.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Done!")
