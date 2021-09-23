mergedTable <- read.table("merge_genes_cobi_and_tnf_annotated.tab", sep="\t", header = TRUE)
head(mergedTable)

# Set fourth column to be the rownames
rownames(mergedTable) <- mergedTable[,1]
head(mergedTable)

mergedTable = as.data.frame(mergedTable)

# plot the log2fc of both treatments
ggplot(mergedTable, aes(log2FoldChange.tnf, log2FoldChange.cobi), colour="grey") +
  scale_color_discrete(name = 'Labels') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Log2FC Cobi vs TNF") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="Log2 Fold Change TNF-alpha", y = "Log2 Fold Change Cobimetinib") +
  # If the genes are not affected by any treatments, colour them grey
  geom_point(data=mergedTable, alpha=0.2, color="grey") +
  # If the genes are affected by Cobi (padj.x<0.05 and log2FoldChange.x<0), change dot color to green
  geom_point(data=mergedTable[which(mergedTable$log2FoldChange.cobi<0 & mergedTable$padj.cobi<0.05),], color="seagreen", alpha=0.2) +
  # If the genes are affected by TNF (padj.y<0.05 and log2FoldChange.y>0), change dot color to red
  geom_point(data=mergedTable[which(mergedTable$log2FoldChange.tnf>0 & mergedTable$padj.tnf<0.05),], color="red3", alpha=0.2) +
  # If the genes are affected by both treatments, colour them blue
  geom_point(data=mergedTable[which(mergedTable$log2FoldChange.tnf>0 & mergedTable$padj.tnf<0.05 & mergedTable$log2FoldChange.cobi<0 & mergedTable$padj.cobi<0.05 & mergedTable$bubbleSize==1),], aes(size=bubbleSize), alpha=0.25, color="midnightblue") +
  geom_point(data=mergedTable[which(mergedTable$log2FoldChange.tnf>0 & mergedTable$padj.tnf<0.05 & mergedTable$log2FoldChange.cobi<0 & mergedTable$padj.cobi<0.05 & mergedTable$bubbleSize==3),], aes(size=bubbleSize), alpha=0.99, fill="midnightblue", pch=21, colour="grey") +
  # Add x and y limits so that 0,0 is in the center
  #xlim(-5, 10) + ylim(-7,10) +
  scale_x_continuous(limits=c(-5, 5), expand = c(0, 0)) +
  scale_y_continuous(limits=c(-5, 5), expand = c(0, 0)) +
  # Add a thicker dashed line at 0
  #geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = 'dashed') 
  geom_segment(aes(x=0,xend=5,y=0,yend=0), size = 0.4, linetype = 'dashed') +
  geom_segment(aes(y=0,yend=-5,x=0,xend=0), size = 0.4, linetype = 'dashed') +
  scale_size_continuous(range = c(2, 4)) +
  #annotate("segment", x = 0.378899136063864, xend = 4, y = -0.419542583178495, yend = -4, colour = "black") +
  #annotate("text", x = 4, y = -4, label = "ATG12") +
  theme(legend.position = "none")
