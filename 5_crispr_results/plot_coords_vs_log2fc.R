# Clear workspace
rm(list = ls())

library(ggplot2)
library(ggrepel)

# Read in the deseq2 table
deseq_table_atg12ltr10i <- read.table("deseq2_results_genes_atg12ltr10i_vs_ctrl.tab", header = TRUE, sep="\t")
head(deseq_table_atg12ltr10i)

# Read in the gene coordinates table and make the gene names the row names
gencode34_genename_genecoords <- read.table("gencode34_genename_genecoords.txt", header = FALSE, sep="\t", row.names=1)
head(gencode34_genename_genecoords)  
nrow(gencode34_genename_genecoords)

# Merge tables by row names to add gene coordinates to each gene name
deseq_with_gene_coords <- merge(deseq_table_atg12ltr10i, gencode34_genename_genecoords, by=0)
head(deseq_with_gene_coords)
nrow(deseq_with_gene_coords)
#14772

# Set first column (gene names) to be the rownames again
rownames(deseq_with_gene_coords) <- deseq_with_gene_coords[,1]
head(deseq_with_gene_coords)

# Add better names to header
colnames(deseq_with_gene_coords) = c("geneName", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "geneChr", "geneStart", "geneEnd") 
head(deseq_with_gene_coords)

# Plot it
ggplot(deseq_with_gene_coords, aes(geneStart, log2FoldChange), colour="grey") +
  scale_color_discrete(name = 'Labels') +
  theme_bw() + 
  geom_hline(aes(yintercept = 0), size=0.2) +
  ggtitle("CRISPRi: LTR10.ATG12 vs Control") +
  theme(plot.title = element_text(hjust = 0.5)) +
  # change plot margins (order is top, right, bottom, left)
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(y="Log2 Fold Change", x = "Genome Coordinates") +
  # Set all dots color to grey
  geom_point(data=deseq_with_gene_coords[which(deseq_with_gene_coords$geneChr=="chr5"),], colour = "grey", size = 1.5) + 
  # add annotation for the targeted LTR10 element
  annotate("rect", xmin=115928578-300000, xmax=115928578+300000, ymin=-0.09, ymax=0.09, alpha=0.99, fill="black") +
  annotate("text", x=115928578, y=0.3, label = "LTR10.ATG12", size=4.5, family="Helvetica") +
  # add annotation for the scale bar
  annotate("rect", xmin=115928578-5000000, xmax=115928578-5000000+1500000, ymin=-3, ymax=-3.02, alpha=0.99, fill="black") +
  annotate("text", x=115928578-5000000+750000, y=-3.2, label = "1.5Mb", size=4.5, family="Helvetica") +
  # If pvalue<0.05, change dot color to blue
  geom_point(data=deseq_with_gene_coords[which(deseq_with_gene_coords$geneChr=="chr5" & deseq_with_gene_coords $padj <0.05 & deseq_with_gene_coords $log2FoldChange>0),],  colour = "blue", size = 3) +
  # Highlight in red genes within 1Mb of targeted LTR10 (still requiring pval<0.05 and log2fc<0)
  geom_point(data=deseq_with_gene_coords[which(deseq_with_gene_coords$geneChr=="chr5" & deseq_with_gene_coords $padj<0.05 & deseq_with_gene_coords $log2FoldChange<0),], colour = "red", size=3) +
  # Add text label for down-regulated genes within 1.5Mb of LTR10 
  geom_text_repel(data=deseq_with_gene_coords[which(deseq_with_gene_coords$geneChr=="chr5" & deseq_with_gene_coords$padj<0.05 & deseq_with_gene_coords$log2FoldChange<0 & deseq_with_gene_coords $geneStart>(115928578-1500000) & deseq_with_gene_coords $geneStart<(115928578+1500000)),], mapping = aes(geneStart, log2FoldChange, label = rownames(deseq_with_gene_coords[which(deseq_with_gene_coords$geneChr=="chr5" & deseq_with_gene_coords $padj<0.05 & deseq_with_gene_coords $log2FoldChange<0 & deseq_with_gene_coords $geneStart>(115928578-1500000) & deseq_with_gene_coords $geneStart<(115928578+1500000)),])),size=4,force=2) +
  # Set the viewing window to be 10Mb centered around LTR10
  scale_x_continuous(labels = scales::comma, limits = c(115928578-5000000,115928578+5000000)) +
  scale_y_continuous(expand = c(0,0),  limits = c(-3.5,1.5)) +
  theme(text=element_text(size=12,family="Helvetica"))

# export pdf 4.5 by 6 inches

# Coordinates of enhancers and genes:
# LTR10.ATG12 enhancer at chr5 115928578
# LTR10.XRCC4 enhancer at chr5 82969145
# LTR10.MEF2D enhancer at chr1 156398365
# LTR10.FGF2 enhancer at chr4 122786289
# LTR10.MCPH1 enhancer at chr8 6648824
# LTR10.KDM6A enhancer at chrX 45005875

# ATG12 gene at chr5 115828199
# FOSL1 gene at chr11 65892048
# JUN gene at chr1 58780790

