#!/bin/bash

## Script for running deeptools 
##
## Example usage:
## bwDir=. bedDir=gencode_tss bedFile=gencode.v34.genes.tss.bed outDir=. sbatch deeptools_heatmap_from_gencode_bed.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=1-00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J deeptools_gencode
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# load modules
module load singularity

# define constant variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img

# set bigwigs and sample labels
wigs=$(ls $bwDir/*.bw | tr "\n" " ")
samplesLabel=$(ls $bwDir/*.bw | xargs -n 1 basename | sed 's/_treat_pileup.bw//g' | tr "\n" " ")

# run the thing
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
--scoreFileName ${wigs} \
--regionsFileName ${bedDir}/${bedFile} \
--beforeRegionStartLength ${windowLeft} \
--afterRegionStartLength ${windowRight} \
--binSize ${binSize} \
--missingDataAsZero \
-o ${outDir}/"BEDFILE"_${bedFile}.mat.gz \
-p ${numCPU}

echo $(date +"[%b %d %H:%M:%S] Plot heatmap...")

# Use "plotHeatmap" to create a png or pdf
zMin=0
yMin=0

singularity exec --bind /scratch/Users $deeptools \
plotHeatmap \
-m $outDir/"BEDFILE"_${bedFile}.mat.gz \
--outFileName $outDir/"BEDFILE"_${bedFile}.png \
--outFileSortedRegions $outDir/"BEDFILE"_${bedFile}.dt.bed \
--outFileNameMatrix $outDir/"BEDFILE"_${bedFile}.matrix.tab \
--sortRegions descend \
--colorMap Blues \
--zMin $zMin --yMin $yMin \
--samplesLabel $samplesLabel \
--regionsLabel $bedFile

echo $(date +"[%b %d %H:%M:%S] Done!")
