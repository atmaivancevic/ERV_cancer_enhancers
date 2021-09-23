library(ggplot2)
library(ggrepel)

deseq2_results_tetranscripts_cobi_vs_ctrl <- read.table("deseq2_results_tetranscripts_cobi_vs_ctrl.tab", sep="\t", header=TRUE)
head(deseq2_results_tetranscripts_cobi_vs_ctrl)

deseq2_results_tetranscripts_cobi_vs_ctrl = as.data.frame(deseq2_results_tetranscripts_cobi_vs_ctrl)

maplot <- ggplot(deseq2_results_tetranscripts_cobi_vs_ctrl, aes(baseMean, log2FoldChange), colour="grey") +
  scale_color_discrete(name = 'Labels') +
  geom_hline(yintercept = 0, size = 0.5, color = 'grey') +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'grey')) +
  theme(axis.text = element_text(size = 14)) +
  labs(y="Log2 Fold Change", x = "Mean of Normalized Family-Level Transcript Counts") +
  # Set all dots color to grey
  geom_point(data=deseq2_results_tetranscripts_cobi_vs_ctrl, colour = "grey", size = 2) + 
  # If pvalue<0.05, change dot color to green
  geom_point(data=deseq2_results_tetranscripts_cobi_vs_ctrl[which(deseq2_results_tetranscripts_cobi_vs_ctrl$padj<0.05),], colour = "red", size = 2) +
  # Add text label for sig points
  geom_text_repel(data = deseq2_results_tetranscripts_cobi_vs_ctrl[which(deseq2_results_tetranscripts_cobi_vs_ctrl$padj<0.05 & deseq2_results_tetranscripts_cobi_vs_ctrl$log2FoldChange<(-1.3)),], mapping = aes(baseMean, log2FoldChange, label = rownames(deseq2_results_tetranscripts_cobi_vs_ctrl[which(deseq2_results_tetranscripts_cobi_vs_ctrl$padj<0.05 & deseq2_results_tetranscripts_cobi_vs_ctrl$log2FoldChange<(-1.3) ),])),size = 5,force = 1) +
  ylim(-2.5, 2.5) + 
  theme(text=element_text(size=14,  family="Helvetica")) 

maplot + scale_x_continuous(trans='log10')
