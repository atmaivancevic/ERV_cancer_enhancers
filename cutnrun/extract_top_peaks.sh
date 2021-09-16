#!/bin/bash

for i in *.bed; do cat $i | sort -nrk5 | head -n 20000 > ${i%.bed}_top20k.bed; done

module load bedtools

cat *_top20k.bed | bedtools sort -i - | bedtools merge -i - -d 100 | awk '{print $0 "\t" "region"NR}' > ALL_CnR_top20k_peaks_merged.bed
