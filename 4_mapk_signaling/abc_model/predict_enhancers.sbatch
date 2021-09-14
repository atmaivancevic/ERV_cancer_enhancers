#!/bin/bash

## Example usage:
## sbatch predict_enhancers.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --time=12:00:00
#SBATCH --mem=64G

# Job name and output
#SBATCH -J predict_enhancers
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

echo $(date +"[%b %d %H:%M:%S] Run predict.py...")

python ABC-Enhancer-Gene-Prediction/src/predict.py \
--enhancers neighborhoods/EnhancerList.txt \
--genes neighborhoods/GeneList.txt \
--HiCdir averageHiC_hg19/ \
--chrom_sizes ABC-Enhancer-Gene-Prediction/reference/chr_sizes \
--hic_resolution 5000 \
--scale_hic_using_powerlaw \
--threshold .02 \
--cellType HCT116 \
--outdir predictions/ \
--make_all_putative

echo $(date +"[%b %d %H:%M:%S] Done!")
