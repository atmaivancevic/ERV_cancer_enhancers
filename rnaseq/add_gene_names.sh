#!/bin/bash

# Example commands to add gene names to feature counts table

# Remove header lines from table
cat featureCounts.txt \
| grep -v "#" \
| grep -v "Geneid" \
> featureCounts_noHeader.txt

# add gene names by matching gene ids
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' featureCounts_noHeader.txt gencode34_geneid_genename.txt \
| paste - featureCounts_noHeader.txt \
> featureCounts_withGeneNames.txt 
