#!/bin/bash

## Script for running deeptools 
##
## Example usage:
## outDir=. outFile=LTR10AF_CnR_H3K27ac sbatch deeptools_H3K27ac_CnR_over_LTR10AF.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=1-00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J deeptools_heatmap
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# load modules
module load singularity

# define constant variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img

# set sample labels
samplesLabel=$(echo Cobi_K27ac_rep1 Cobi_K27ac_rep2 UT_K27ac_rep1 UT_K27ac_rep2 TNFa_K27ac_rep1 TNFa_K27ac_rep2 | tr " " "\n")

# run the thing
pwd; hostname; date

echo "Starting deeptools..."
echo $(date +"[%b %d %H:%M:%S] Compute matrix...")

# Use "computeMatrix" to generate data underlying heatmap
windowLeft=4000
windowRight=4000
binSize=10
numCPU=8

singularity exec --bind /Shares/CL_Shared $deeptools \
computeMatrix reference-point \
--referencePoint TSS \
--scoreFileName Cobi_K27ac_rep1_treat_pileup.bw \
Cobi_K27ac_rep2_treat_pileup.bw \
UT_K27ac_rep1_treat_pileup.bw \
UT_K27ac_rep2_treat_pileup.bw \
TNF_K27ac_rep1_treat_pileup.bw \
TNF_K27ac_rep2_treat_pileup.bw \
--regionsFileName LTR10A_and_LTR10F_merged2kb.bed \
--beforeRegionStartLength ${windowLeft} \
--afterRegionStartLength ${windowRight} \
--binSize ${binSize} \
--missingDataAsZero \
-o ${outDir}/$outFile.mat.gz \
-p ${numCPU}

echo $(date +"[%b %d %H:%M:%S] Plot heatmap...")

# Use "plotHeatmap" to create a png or pdf
zMin=0
yMin=0

singularity exec --bind /Shares/CL_Shared $deeptools \
plotHeatmap \
-m $outDir/$outFile.mat.gz \
--outFileName $outDir/$outFile.pdf \
--outFileSortedRegions $outDir/$outFile.dt.bed \
--outFileNameMatrix $outDir/$outFile.matrix.tab \
--sortRegions descend \
--colorMap Blues \
--zMin $zMin --yMin $yMin \
--samplesLabel $samplesLabel \
--regionsLabel LTR10AF_merged

echo $(date +"[%b %d %H:%M:%S] Done!")
