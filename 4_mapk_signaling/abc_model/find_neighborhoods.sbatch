#!/bin/bash

## Example usage:
## sbatch find_neighborhoods.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --time=12:00:00
#SBATCH --mem=64G

# Job name and output
#SBATCH -J run_neighborhoods
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Run the thing
pwd; hostname; date

module purge
module load gcc/7.2.0
module load bedtools

echo $(date +"[%b %d %H:%M:%S] Initialize the shell to use conda...")
eval "$(conda shell.bash hook)"

echo $(date +"[%b %d %H:%M:%S] Activate the ABC conda environment...")
conda activate final-abc-env

python ABC-Enhancer-Gene-Prediction/src/run.neighborhoods.py \
--candidate_enhancer_regions HCT116_ATACseq_SRR8544480_peaks.narrowPeak.sorted.candidateRegions.bed \
--genes ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.bed \
--H3K27ac Untreated_K27ac.sorted.bam \
--ATAC HCT116_ATACseq_SRR8544480_trimmed.sorted.uniq.bam \
--chrom_sizes ABC-Enhancer-Gene-Prediction/reference/chr_sizes \
--ubiquitously_expressed_genes ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19.txt \
--cellType HCT116 \
--outdir neighborhoods

echo $(date +"[%b %d %H:%M:%S] Done!")
