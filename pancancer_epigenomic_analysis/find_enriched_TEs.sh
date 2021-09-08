#!/bin/bash

# concatenate all healthy adult regulatory regions from Roadmap
module load bedtools

cat *1_TssA.bed *6_EnhG.bed *7_Enh.bed \
| bedtools sort -i - \
| bedtools merge -i - -c 4 -o distinct \
> healthy_adult_regions.bed

# for each cancer, subtract the healthy regulatory regions
for i in *_peakCalls.bed; 
	do bedtools subtract -a $i -b healthy_adult_regions.bed -A > ${i%.bed}_minus_healthy.bed; 
	done

# for each cancer, bgzip the output bed file (required for GIGGLE)
module load samtools

for i in *peakCalls_minus_healthy.bed; 
	do bgzip $i; 
	done

# for each cancer, search the predicted cancer-specific regulatory regions against human TE families
for i in *peakCalls_minus_healthy.bed.gz; 
	do echo $i; giggle search -q $i -i giggle/hg38/repeats/indexed -s -g 3209286105 \
	| sed 's#sorted/##g' \
	| sed 's/.bed.gz//g' \
	| grep -v "#" \
	| sort -nrk8 \
	| sed '1i#file\tfile_size\toverlaps\todds_ratio\tfishers_two_tail\tfishers_left_tail\tfishers_right_tail\tcombo_score' \
	> ${i%.gz}.VSrepeats.tab; 
	done

# combine all cancers and all TEs into one table 
# and rank by highest to lowest GIGGLE enrichment score

for i in *_minus_healthy.bed.VSrepeats.tab; 
	do cat $i \
	| grep -v "#" \
	| awk '{print "'$i'" "\t" $0}' \
	| sed 's/_peakCalls_minus_healthy.bed.VSrepeats.tab//g' \
	| sed '1i#cancer\trepeat\tfile_size\toverlaps\todds_ratio\tfishers_two_tail\tfishers_left_tail\tfishers_right_tail\tcombo_score' \
	> $i.tmp; 
	done

cat *.tmp \
| grep -v "#" \
| sort -nrk9 \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $9  }' \
| sed '1icancer\trepeat\tfilesize\toverlaps\toddsratio\tgigglescore' \
> tcgaMinusRoadmap_vs_TEs_rankedByScore.tab

##################################################################################

# The following steps are for visualization purposes
### 1. make a pan-cancer bubble plot showing the top TE-cancer associations

# remove satellite and low complexity repeats
cat tcgaMinusRoadmap_vs_TEs_rankedByScore.tab \
| grep -v "G-rich" \
| grep -v "GA-rich" \
| grep -v "U1" \
| grep -v "BSR_Beta" \
| grep -v "GSATII" \
| grep -v "GSATX" \
| grep -v "MSR1" \
| grep -v "REP522" \
| grep -v "TAR1" \
> tcgaMinusRoadmap_vs_TEs_rankedByScore_noLowComplxRepeats.tab

# set significance thresholds for TE-cancer associations
# i.e. must have over 25 overlaps, odds ratio over 10, and giggle score over 100 in at least one cancer
cat tcgaMinusRoadmap_vs_TEs_rankedByScore_noLowComplxRepeats.tab \
| awk '{if ($4>25) print}' \
| awk '{if ($5>10) print}' \
| awk '{if ($6>100) print}' \
> tcgaMinusRoadmap_vs_TEs_rankedByScore_noLowComplxRepeats_overlaps25_oddsratio10_score100.tab

# this leaves 23 TE families
# extract the giggle results for these 23 top TEs
cat tcgaMinusRoadmap_vs_TEs_rankedByScore_noLowComplxRepeats.tab \
| awk '{if ($2=="HERVE_a-int"|| $2=="Harlequin-int"|| $2=="LTR10A"|| $2=="LTR10C"|| $2=="LTR10F"|| $2=="LTR13"|| $2=="LTR13A"|| $2=="LTR2B"|| $2=="LTR3A"|| $2=="LTR3B"|| $2=="LTR5B"|| $2=="LTR5_Hs"|| $2=="LTR7Y"|| $2=="LTR9"|| $2=="MER11A"|| $2=="MER11B"|| $2=="MER11D"|| $2=="MER39B"|| $2=="MER41B"|| $2=="MER44B"|| $2=="MER51A"|| $2=="PABL_A"|| $2=="SVA_A") print}' \
| sed 1i"cancer_type\trepeat_family\trepeat_family_size\toverlaps\todds_ratio\tgiggle_score" \
> SuppTable1_TCGA_giggle_results_top23TEs.tab

# import and plot in jupyter (Fig1_bubbles.py)

### 2. make volcano plots showing top enriched TEs for each cancer type
# e.g. prepare the colon cancer (COAD) table for plotting in R

cat COAD_peakCalls_minus_healthy.bed.VSrepeats.tab \
| sed 's/#//g' \
| awk '{print $1 "\t" $4 "\t" $5 "\t" $8}' \
| sed 's/file/repeat/g' \
| sed 's/fishers_two_tail/fishers_two_tail_pval/g' \
> COAD_minus_Roadmap.tab.volcano

# import and plot in R (Fig1_volcano_COAD.R)
