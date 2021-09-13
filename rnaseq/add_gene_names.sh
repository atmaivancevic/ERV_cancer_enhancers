#!/bin/bash

# Example commands to add gene names to feature counts table

# Remove header lines from table
cat featureCounts.txt \
| grep -v "#" \
| grep -v "Geneid" \
> featureCounts_noHeader.txt

# Add gene names by matching gene ids
# And add header back in
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' featureCounts_noHeader.txt gencode34_geneid_genename.txt \
| paste - featureCounts_noHeader.txt \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15}' \
| sed '1iGeneid\tGeneName\tGeneDesc\tLength\tCobi_24_1_S54\tCobi_24_2_S55\tTNF_24_1_S58\tTNF_24_2_S59\tUT_24_1_S52\tUT_24_2_S53' \
> featureCounts_withGeneNames.txt 
