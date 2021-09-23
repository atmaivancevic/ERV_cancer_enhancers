norm_counts_74_genes <- read.table("74_genes_normalized_counts_rnaseq.tab", sep="\t", header = TRUE)
head(norm_counts_74_genes)

# Set first column to be the rownames
rownames(norm_counts_74_genes) <- norm_counts_74_genes[,1]
head(norm_counts_74_genes)

# Remove the first column
norm_counts_74_genes <- norm_counts_74_genes[ ,2:ncol(norm_counts_74_genes)]
head(norm_counts_74_genes)

norm_counts_74_genes <- norm_counts_74_genes[ ,c("Cobi_24h_1", "Cobi_24h_2", "UT_24h_1", "UT_24h_2", "TNF_24h_1", "TNF_24h_2")]
head(norm_counts_74_genes)

norm_counts_74_genes = as.data.frame(norm_counts_74_genes)

### Run pheatmap
pheatmap(norm_counts_74_genes, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         cluster_rows = T, clustering_distance_rows = "manhattan", show_rownames=T, fontsize = 10,
         fontsize_row = 10, height=20, cluster_cols = F, scale="row")
