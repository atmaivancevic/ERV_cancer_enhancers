#!/bin/bash

## Script for running featureCounts
##
## Example usage:
## inDir=. outDir=. gtfFile=/Shares/CL_Shared/db/genomes/hg38/annotations/gencode.v34.annotation.gtf feature=exon strandOption=2 sbatch feature_counts.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=2:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J featureCounts
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# load modules
module load subread

# Run featureCounts, outputting a single table containing count information for all samples
pwd; hostname; date

echo "Starting featureCounts..."
echo "featureCounts version: "$(featureCounts -v)

inputBams=${inDir}/*.uniq.bam

featureCounts -T 1 -O -p -s $strandOption -t $feature -g gene_id -a ${gtfFile} -o ${outDir}/featureCounts.txt $inputBams

echo $(date +"[%b %d %H:%M:%S] Done!")

## Explanation of arguments:
# '-O' - assign reads to all their overlapping meta-features
# '-T <int>' - number of threads; 1 by default
# '-s <int>' - perform strand-specific read counting; options: 0 (unstranded), 1 (stranded), and 2 (reversely stranded); 0 by default
# '-g <option>' - specify attribute type in GTF annotation; 'gene_id' by default
# '-a <file.gtf' - name of annotation file; gtf format by default
