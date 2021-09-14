#!/bin/bash

module load bedtools
cd predictions/

# Convert bedpe file to bed
gunzip EnhancerPredictionsAllPutative.txt.gz
cat EnhancerPredictionsAllPutative.txt \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $7 "\t" $21 "\t" "."}' \
| grep -v ABC.Score \
| bedtools sort -i - \
> EnhancerPredictionsAllPutative.sorted.bed

# Columns are now: 
# enhancerChr    enhancerStart    enhancerStop    geneName   ABCscore    strand

cd ..
mkdir predictions_lifted_to_hg38

liftOver predictions/EnhancerPredictionsAllPutative.sorted.bed hg19ToHg38.over.chain predictions_lifted_to_hg38/EnhancerPredictionsAllPutative.hg38.bed predictions_lifted_to_hg38/EnhancerPredictionsAllPutative.hg19Tohg38Unmapped.bed

cd predictions_lifted_to_hg38/
# 10692 EnhancerPredictionsAllPutative.hg19Tohg38Unmapped.bed
# 9765253 EnhancerPredictionsAllPutative.hg38.bed
