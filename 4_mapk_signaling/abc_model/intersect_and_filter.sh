#!/bin/bash

module load bedtools

# Intersect predicted enhancers with LTR10A/F elements
bedtools intersect -a EnhancerPredictionsAllPutative.hg38.bed -b LTR10_merged.bed -loj \
| grep LTR10 \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "." "\t" $10}' \
| uniq \
| grep 'LTR10A\|LTR10F' \
> EnhancerPredictionsAllPutativehg38_intersect_LTR10AF.bed
# 8780

# Filter for enhancer-gene pairs with ABC interaction score over 0.001
cat EnhancerPredictionsAllPutativehg38_intersect_LTR10AF.bed \
| awk '{if ($5!="NaN") print}' \
| awk '{if ($5>0.001) print}' \
| sort -nrk5 \
> EnhancerPredictionsAllPutativehg38_intersect_LTR10AF_scoreOver0.001.bed
# 3583
