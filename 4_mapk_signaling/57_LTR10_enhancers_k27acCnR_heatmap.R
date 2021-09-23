# Clear workspace
rm(list = ls())

# Read table
norm_counts_57_LTR10s <- read.table("57_LTR10_enhancers_normalized_counts_h3K27ac_cutnrun.tab", sep="\t", header = TRUE)
head(norm_counts_57_LTR10s)

# Set first column to be the rownames
rownames(norm_counts_57_LTR10s) <- norm_counts_57_LTR10s[,1]
head(norm_counts_57_LTR10s)

# Remove the first column
norm_counts_57_LTR10s <- norm_counts_57_LTR10s[ ,2:ncol(norm_counts_57_LTR10s)]
head(norm_counts_57_LTR10s)

# Reorder the columns
norm_counts_57_LTR10s <- norm_counts_57_LTR10s[, c("Cobi_K27ac_rep1", "Cobi_K27ac_rep2", "UT_K27ac_rep1", "UT_K27ac_rep2", "TNFa_K27ac_rep1", "TNFa_K27ac_rep2")]
head(norm_counts_57_LTR10s)

norm_counts_57_LTR10s = as.data.frame(norm_counts_57_LTR10s)

### Plot a heatmap
library(RColorBrewer)
library(pheatmap)
pheatmap(norm_counts_57_LTR10s, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         cluster_rows = T, cluster_cols = F, show_rownames=T, fontsize = 10, clustering_distance_rows = "manhattan",
         fontsize_row = 10, height=20, scale="row")
