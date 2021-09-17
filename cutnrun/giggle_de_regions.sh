#!/bin/bash

# Extract significantly DE regions for each comparison

# Cobi vs untreated

cat deseq2_results_H3K27ac_CnR_cobi_vs_ctrl.tab \
| grep -v baseMean \
| awk '{if ($6<0.05 && $3<0) print}' \
| awk '{print $1}' \
| sed 's/:/\t/g' \
| sed 's/-/\t/g' \
> DE_regions_H3K27ac_CnR_cobi_vs_ctrl.bed

# Tnf vs untreated

cat deseq2_results_H3K27ac_CnR_tnf_vs_ctrl.tab \
| grep -v baseMean \
| awk '{if ($6<0.05 && $3>0) print}' \
| awk '{print $1}' \
| sed 's/:/\t/g' \
| sed 's/-/\t/g' \
> DE_regions_H3K27ac_CnR_tnf_vs_ctrl.bed

# Bgzip and giggle search

module load samtools

for i in *_ctrl.bed; do bgzip $i; done

giggle search -q DE_regions_H3K27ac_CnR_cobi_vs_ctrl.bed.gz \
-i /Shares/CL_Shared/db/giggle/hg38/repeats/indexed \
-s -g 3209286105 \
| grep -v "#" \
| sort -nrk8 \
| sed 's#sorted/##g' \
| sed 's/.bed.gz//g' \
| sed '1irepeat_family\trepeat_family_size\toverlaps\todds_ratio\tfishers_two_tail\tfishers_left_tail\tfishers_right_tail\tgiggle_score' \
> DE_regions_H3K27ac_CnR_cobimetinib_giggle_results.tab

giggle search -q DE_regions_H3K27ac_CnR_tnf_vs_ctrl.bed.gz \
-i /Shares/CL_Shared/db/giggle/hg38/repeats/indexed \
-s -g 3209286105 \
| grep -v "#" \
| sort -nrk8 \
| sed 's#sorted/##g' \
| sed 's/.bed.gz//g' \
| sed '1irepeat_family\trepeat_family_size\toverlaps\todds_ratio\tfishers_two_tail\tfishers_left_tail\tfishers_right_tail\tgiggle_score' \
> DE_regions_H3K27ac_CnR_tnfalpha_giggle_results.tab
