# ABC output file:
# 1) enhancer predictions: EnhancerPredictionsAllPutativehg38_intersect_LTR10AF_scoreOver0.001.bed

# MAPK output files:
# 1) Merged table of Cobi/TNF-alpha genes: merge_genes_cobi_and_tnf.tab
# 2) Normalized counts from Cut&Run H3K27ac regions: normalized_counts_H3K27ac_CnR.tab

module load bedtools

# How many unique LTR10-derived predicted enhancers? 
cat EnhancerPredictionsAllPutativehg38_intersect_LTR10AF_scoreOver0.001.bed \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $7}' \
| bedtools sort -i - \
| bedtools merge -i - -c 4 -o distinct \
> EnhancerPredictionsAllPutativehg38_intersect_LTR10AF_scoreOver0.001_unique.bed

# 134 distinct LTR10s 

# How many of these have a H3K27ac peak in HCT116 cells?
# Turn the normalized count table of H3K27ac regions into a bed file
# Then intersect with predicted enhancers

cat normalized_counts_H3K27ac_CnR.tab | grep -v UT_K27ac_rep1 | awk '{if ($2+$3+$4+$5+$6+$7!=0) print}' | awk '{print $1}' | sed 's/:/\t/g' | sed 's/-/\t/g' > H3K27ac_CnR_nonzero_regions.bed

bedtools intersect -a EnhancerPredictionsAllPutativehg38_intersect_LTR10AF_scoreOver0.001_unique.bed -b H3K27ac_CnR_nonzero_regions.bed -loj | grep -v "\-1" | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > 98_LTR10_derived_enhancers.bed 

# 98 LTR10s left

# How many of these 98 LTR10s are near MAPK-regulated genes?
# e.g. within 1.5Mb of a gene that showed a significant change with both MAPK treatments 
# ([Cobi_Log2FC<0 & Cobi_padj<0.05] AND [TNFLog2FC>0 & TNF_padj<0.05])?

# add gene coordinates to the gene names in the merged tnf/cobi deseq2 table
awk 'FNR==NR{a[$1]=$1"\t"$2"\t"$3"\t"$4"\t"$5;next}{print $0 "\t" a[$1]}' gencode34_genename_genecoords.txt  merge_genes_cobi_and_tnf.tab | grep -v log2FoldChange.x | awk '{print $14 "\t" $15 "\t" $16 "\t" $13 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}' | sed '1ichr\tstart\tstop\tgene\tbaseMean.Cobi\tlog2FoldChange.Cobi\tlfcSE.Cobi\tpvalue.Cobi\tpadj.Cobi\tbaseMean.tnf\tlog2FoldChange.tnf\tlfcSE.tnf\tpvalue.tnf\tpadj.tnf' > merge_genes_cobi_and_tnf_with_genecoords.tab

# then look for LTR10 enhancers within 1.5Mb of MAPK-regulated genes
cat merge_genes_cobi_and_tnf_with_genecoords.tab | grep -v baseMean.Cobi > merge_genes_cobi_and_tnf_with_genecoords_noheader.bed

bedtools window -a 98_LTR10_derived_enhancers.bed -b merge_genes_cobi_and_tnf_with_genecoords_noheader.bed -w 1500000 | awk '{if (($10<0 && $13<0.05) && ($15>0 && $18<0.05)) print}' | awk '{print $8}' | sort | uniq | wc -l
# 74 distinct genes

bedtools window -a 98_LTR10_derived_enhancers.bed -b merge_genes_cobi_and_tnf_with_genecoords_noheader.bed -w 1500000 | awk '{if (($10<0 && $13<0.05) && ($15>0 && $18<0.05)) print}' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' | sort | uniq | wc -l
# 57 distinct LTR10 enhancers

# make a candidate list containing these final enhancer-gene pairs
bedtools window -a 98_LTR10_derived_enhancers.bed -b merge_genes_cobi_and_tnf_with_genecoords_noheader.bed -w 1500000 | awk '{if (($10<0 && $13<0.05) && ($15>0 && $18<0.05)) print}' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $8}' | bedtools merge -i - -c 5 -o distinct > LTR10_predicted_enhancer_gene_pairs.bed
