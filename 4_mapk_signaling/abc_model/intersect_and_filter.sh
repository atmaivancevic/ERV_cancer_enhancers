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

########## UP TO HERE! 


# Intersect with H3K27ac-marked LTR10A/F only
bedtools intersect -a EnhancerPredictionsAllPutativehg38_intersect_LTR10AF_scoreOver0.001_unique.bed -b normalized_counts_H3K27ac_CnR_modified.tmp -loj | grep -v "\-1"  | awk '{if ($8!="0") print}' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > 98_LTR10_derived_enhancers.bed

