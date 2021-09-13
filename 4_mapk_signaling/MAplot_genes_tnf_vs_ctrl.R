# MA plot

# Read table
deseq2_results_tnf_vs_ctrl <- read.table("deseq2_results_tnf_vs_ctrl.tab", sep="\t", header = TRUE)
head(deseq2_results_tnf_vs_ctrl)

# Convert to data frame
deseq2_results_tnf_vs_ctrl = as.data.frame(deseq2_results_tnf_vs_ctrl)

# Plot
plot1 <- ggplot(deseq2_results_tnf_vs_ctrl, aes(baseMean, log2FoldChange), colour="grey") +
  scale_color_discrete(name = 'Labels') +
  theme_bw() +
  ggtitle("MA plot: tnf_24hr vs ctrl_24hr") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y="Log2 Fold Change", x = "Mean of Normalized Counts") +
  # Set all dots color to grey
  geom_point(data=deseq2_results_tnf_vs_ctrl, colour = "grey", alpha=0.5) +
  # If padj<0.05 and log2fc<0, change dot color to red
  geom_point(data=deseq2_results_tnf_vs_ctrl[which(deseq2_results_tnf_vs_ctrl$log2FoldChange>0 & deseq2_results_tnf_vs_ctrl$padj<0.05),], colour = "red3", , alpha=0.5) +
  # Add text labels for some points
  geom_text_repel(data = deseq2_results_tnf_vs_ctrl[which(deseq2_results_tnf_vs_ctrl$padj<0.00000000000000000000000000000000000000000000000000000000000000005),], mapping = aes(baseMean, log2FoldChange, label = rownames(deseq2_results_tnf_vs_ctrl[which(deseq2_results_tnf_vs_ctrl$padj<0.00000000000000000000000000000000000000000000000000000000000000005),])),size = 4,force = 1)

plot1 + scale_x_continuous(trans='log10')
