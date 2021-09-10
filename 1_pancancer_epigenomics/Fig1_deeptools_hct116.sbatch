#!/bin/bash

## Script for running deeptools to generate a heatmap
##
## Example usage:
## outDir=. outFile=LTR10AF_over_HCT116_marks sbatch Fig1_deeptools_hct116.sbatch

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

# Load modules
module load singularity

# Define constant variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img

# Set sample labels
samplesLabel=$(echo H3K27ac_84767 H3K4me1_85458 POLR2A_46201 EP300_42907 ATACseq_105341 | tr " " "\n")

# Run the thing
pwd; hostname; date

echo "Starting deeptools..."
echo $(date +"[%b %d %H:%M:%S] Compute matrix...")

# Use "computeMatrix" to generate data underlying heatmap
windowLeft=4000
windowRight=4000
binSize=10
numCPU=8

singularity exec --bind /scratch/Users $deeptools \
computeMatrix reference-point \
--referencePoint TSS \
--scoreFileName HCT116_H3K27ac_84767.bw HCT116_H3K4me1_85458.bw HCT116_POLR2A_46201.bw HCT116_EP300_42907.bw HCT116_ATACseq_105341.bw \
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

singularity exec --bind /scratch/Users $deeptools \
plotHeatmap \
-m $outDir/$outFile.mat.gz \
--outFileName $outDir/$outFile.pdf \
--outFileSortedRegions $outDir/$outFile.dt.bed \
--outFileNameMatrix $outDir/$outFile.matrix.tab \
--sortRegions descend \
--colorMap Blues \
--zMin $zMin --yMin $yMin \
--yMax 0.6 0.6 0.3 0.3 0.3 \
--samplesLabel $samplesLabel \
--regionsLabel LTR10AF_merged

echo $(date +"[%b %d %H:%M:%S] Done!")
