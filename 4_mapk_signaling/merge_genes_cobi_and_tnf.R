# Merge the TNF-alpha and Cobimetinib results

# Read tables
deseq2_results_genes_cobi_vs_ctrl <- read.table("deseq2_results_genes_cobi_vs_ctrl.tab", sep="\t", header = TRUE)
deseq2_results_genes_tnf_vs_ctrl <- read.table("deseq2_results_genes_tnf_vs_ctrl.tab", sep="\t", header = TRUE)

head(deseq2_results_genes_cobi_vs_ctrl)
head(deseq2_results_genes_tnf_vs_ctrl)

# Merge the tables by common gene names
# by=0 means merge by row names
merge_genes_cobi_and_tnf <- merge(deseq2_results_genes_cobi_vs_ctrl, deseq2_results_genes_tnf_vs_ctrl, by=0)
head(merge_genes_cobi_and_tnf)
nrow(merge_genes_cobi_and_tnf)
#14790

rownames(merge_genes_cobi_and_tnf) <- merge_genes_cobi_and_tnf[,1]
head(merge_genes_cobi_and_tnf)

# Sort by lowest padj value in Cobimetinib (cobi is x, TNF is y)
merge_genes_cobi_and_tnf <- merge_genes_cobi_and_tnf[order(merge_genes_cobi_and_tnf$padj.x), ]
head(merge_genes_cobi_and_tnf, 10)

# Which genes show a sig diff with both treatments?
our_genes = subset(merge_genes_cobi_and_tnf, log2FoldChange.y>0 & padj.y<0.05 & log2FoldChange.x<0 & padj.x<0.05)
nrow(our_genes)
our_genes
#620

# Save the merged tnf-cobi table
write.table(merge_genes_cobi_and_tnf, file="merge_genes_cobi_and_tnf.tab", quote = FALSE, row.names = TRUE, sep = "\t")
