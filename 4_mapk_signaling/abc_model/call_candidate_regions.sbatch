#!/bin/bash

## Example usage:
## sbatch call_candidate_regions.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=12:00:00
#SBATCH --mem=64G

# Job name and output
#SBATCH -J make_regions
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

echo $(date +"[%b %d %H:%M:%S] Run makeCandidateRegions.py...")
python ABC-Enhancer-Gene-Prediction/src/makeCandidateRegions.py \
--narrowPeak HCT116_ATACseq_SRR8544480_peaks.narrowPeak.sorted \
--bam HCT116_ATACseq_SRR8544480_trimmed.sorted.uniq.bam \
--outDir candidate_regions_output/ \
--chrom_sizes /Shares/CL_Shared/db/genomes/hg19/fa/hg19.chrom.sizes \
--regions_blocklist ABC-Enhancer-Gene-Prediction/reference/wgEncodeHg19ConsensusSignalArtifactRegions.bed \
--peakExtendFromSummit 250 \
--nStrongestPeaks 150000

echo $(date +"[%b %d %H:%M:%S] Done!")
