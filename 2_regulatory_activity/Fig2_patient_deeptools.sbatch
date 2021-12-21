#!/bin/bash

## Script for running deeptools 
##
## Example usage:
## bwDir=. bedDir=. bedFile=LTR10A_and_LTR10F_merged2kb.bed outDir=. outFile=LTR10AF_H3K27ac_H3K4me1_P52_P58 sbatch Fig2_patient_deeptools.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=1-00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J deeptools2
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# load modules
module load singularity

# define constant variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img

# set bigwigs and sample labels
wigs=$(cat bigwigs2.txt)
samplesLabel=$(echo AKCC52_Tumor_K27ac AKCC52_Norm_K27ac AKCC58_Tumor_K27ac AKCC58_Norm_K27ac AKCC52_Tumor_K4me1 AKCC52_Norm_K4me1 AKCC58_Tumor_K4me1 AKCC58_Norm_K4me1 | tr " " "\n")

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
--scoreFileName ${wigs} \
--regionsFileName ${bedDir}/${bedFile} \
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
--zMax 15 15 15 15 35 35 35 35 \
--yMax 8 8 8 8 25 25 25 25 \
--samplesLabel $samplesLabel \
--regionsLabel $bedFile 
#--sortUsingSamples 1

echo $(date +"[%b %d %H:%M:%S] Done!")
