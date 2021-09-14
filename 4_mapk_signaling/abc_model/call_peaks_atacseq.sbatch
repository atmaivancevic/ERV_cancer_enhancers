#!/bin/bash

## Example usage:
## sbatch call_peaks_atacseq.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=12:00:00
#SBATCH --mem=64G

# Job name and output
#SBATCH -J macs2_call_peaks
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Load modules
module load python/2.7.14/MACS/2.1.1

# Run the thing
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Starting macs2...")

macs2 callpeak \
--treatment HCT116_ATACseq_SRR8544480_trimmed.sorted.uniq.bam \
--name HCT116_ATACseq_SRR8544480 \
--format BAM \
-g hs \
--pvalue .1 \
--call-summits \
--outdir .

echo $(date +"[%b %d %H:%M:%S] Done!")
