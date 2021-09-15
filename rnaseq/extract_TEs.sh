#!/bin/bash

for i in *.cntTable; 
	do cat $i \
	| grep -v ENSG00000 \
	| grep -v ":Satellite" \
	| grep -v ":RNA" \
	| grep -v ":Unknown" \
	| grep -v "RTE-Bo\.B" \
	| grep -v "?" \
	| awk -F ":" '{print $1 "\t" $3}' \
	| sed 's#gene/TE#TE\tTE_family#' \
	| awk '{print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $2}' \
	> ${i%.cntTable}_TEsonly.counttab; 
	done
