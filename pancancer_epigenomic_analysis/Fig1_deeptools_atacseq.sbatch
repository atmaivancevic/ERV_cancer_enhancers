#!/bin/bash

## Script for running deeptools to generate a heatmap
##
## Example usage:
## outDir=. outFile=TCGA_ATACseq_over_LTR10AF sbatch Fig1_deeptools_atacseq.sbatch

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

# Load modules
module load singularity

# Define constant variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img

# Set sample labels
samplesLabel=$(echo COAD_0914606C COAD_0B139EBC COAD_0B7D5055 COAD_1E46B832 COAD_2EB7D36F COAD_32AA8B86 COAD_3ECCFFD0 COAD_42D5A2CF | tr " " "\n")

# Run the thing
pwd; hostname; date

echo "Starting deeptools..."
echo $(date +"[%b %d %H:%M:%S] Compute matrix...")

# Use "computeMatrix" to generate data underlying heatmap
windowLeft=4000
windowRight=4000
binSize=10
numCPU=8

singularity exec --bind /Shares/CL_Shared \
$deeptools computeMatrix reference-point \
--referencePoint TSS \
--scoreFileName COAD_0914606C_2CA1_4287_B530_DB70EA93ED6C_X022_S02_L050_B1_T1_P053.insertions.bw \
COAD_0B139EBC_D372_4EC3_90DB_4CC9BC6F38DC_X006_S05_L038_B1_T1_P012.insertions.bw \
COAD_0B7D5055_6D9C_46BB_A420_8098CD555089_X003_S07_L031_B1_T1_P002.insertions.bw \
COAD_1E46B832_4B41_4CA3_B2DB_C8287EFFC589_X011_S04_L007_B1_T1_P025.insertions.bw \
COAD_2EB7D36F_5616_460C_91C9_6D0AEA893D97_X001_S04_L007_B1_T1_P004.insertions.bw \
COAD_32AA8B86_5414_4552_B946_E2F368618E5F_X009_S12_L047_B1_T1_P016.insertions.bw \
COAD_3ECCFFD0_0915_4151_9B10_6C7C37372A72_X002_S05_L009_B1_T1_P001.insertions.bw \
COAD_42D5A2CF_D448_4072_8278_98B643DAA870_X022_S08_L061_B1_T1_P049.insertions.bw \
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

singularity exec --bind /Shares/CL_Shared \
$deeptools plotHeatmap \
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
