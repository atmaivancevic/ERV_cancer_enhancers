# Compare the TNF-alpha and Cobimetinib results directly
# E.g. make a log2FC_Cobi vs log2FC_TNF graph

# Read tables
deseq2_results_cobi_vs_ctrl <- read.table("deseq2_results_cobi_vs_ctrl.tab", sep="\t", header = TRUE)
deseq2_results_tnf_vs_ctrl <- read.table("deseq2_results_tnf_vs_ctrl.tab", sep="\t", header = TRUE)

head(deseq2_results_cobi_vs_ctrl)
head(deseq2_results_tnf_vs_ctrl)

# Merge the tables by common gene names
# by=0 means merge by row names
merged_tnf_cobi <- merge(deseq2_results_cobi_vs_ctrl, deseq2_results_tnf_vs_ctrl, by=0)
head(merged_tnf_cobi)
nrow(merged_tnf_cobi)
#14790

rownames(merged_tnf_cobi) <- merged_tnf_cobi[,1]
head(merged_tnf_cobi)

# Sort by lowest padj value in Cobimetinib (cobi is x, TNF is y)
merged_tnf_cobi <- merged_tnf_cobi[order(merged_tnf_cobi$padj.x), ]
head(merged_tnf_cobi, 10)

# Which genes show a sig diff with both treatments?
our_genes = subset(merged_tnf_cobi, log2FoldChange.y>0 & padj.y<0.05 & log2FoldChange.x<0 & padj.x<0.05)
nrow(our_genes)
our_genes
#620

# Plot the log2FC against each other
ggplot(merged_tnf_cobi, aes(log2FoldChange.y, log2FoldChange.x), colour="grey") +
  scale_color_discrete(name = 'Labels') +
  theme_bw() +
  ggtitle("Log2FC Cobi vs TNF") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="Log2 Fold Change TNF-alpha", y = "Log2 Fold Change Cobimetinib") +
  # Set all dots color to grey
  geom_point(data=merged_tnf_cobi, colour = "grey") +
  # If the genes are affected by Cobi (padj.x<0.05 and log2FoldChange.x<0), change dot color to green
  geom_point(data=merged_tnf_cobi[which(merged_tnf_cobi$log2FoldChange.x<0 & merged_tnf_cobi$padj.x<0.05),], colour = "seagreen") +
  # If the genes are affected by TNF (padj.y<0.05 and log2FoldChange.y>0), change dot color to blue
  geom_point(data=merged_tnf_cobi[which(merged_tnf_cobi$log2FoldChange.y>0 & merged_tnf_cobi$padj.y<0.05),], colour = "red3") +
  # If the genes are affected by both treatments, colour them red
  geom_point(data=merged_tnf_cobi[which(merged_tnf_cobi$log2FoldChange.y>0 & merged_tnf_cobi$padj.y<0.05 & merged_tnf_cobi$log2FoldChange.x<0 & merged_tnf_cobi$padj.x<0.05),], colour = "midnightblue") +
  # Add x and y limits so that 0,0 is in the center
  #xlim(-5, 10) + ylim(-7,10) +
  scale_x_continuous(limits=c(-5, 10), expand = c(0, 0)) +
  scale_y_continuous(limits=c(-7.5, 10), expand = c(0, 0)) +
  # Add a thicker dashed line at 0
  #geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = 'dashed')
  geom_segment(aes(x=0,xend=10,y=0,yend=0), size = 0.4, linetype = 'dashed') +
  geom_segment(aes(y=0,yend=-7.5,x=0,xend=0), size = 0.4, linetype = 'dashed')

# Save the merged tnf-cobi table
write.table(merged_tnf_cobi, file="merged_tnf_cobi.tab", quote = FALSE, row.names = TRUE, sep = "\t")
