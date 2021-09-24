library("DESeq2")
library("pheatmap")
library("ggplot2")
library("ggrepel")
library("reshape2")

# Read in the output from featureCounts
countdata <- read.table("feature_counts_crispri.tab", header = TRUE, sep="\t")
head(countdata)

# Set second column (gene names) to be the rownames
rownames(countdata) <- countdata[,2]
head(countdata)

# Remove the first two columns (geneid, genename)
countdata <- countdata[ ,3:ncol(countdata)]
head(countdata)

# Convert countdata table into a matrix - necessary for running DESeq2. 
countdata <- as.matrix(countdata)
head(countdata)

# Assign control vs test samples
condition <- factor(c(rep("ATG12_LTR10i", 2), rep("ATG12_TSSi", 2), rep("cJun_TSSi", 2), rep("FGF2_LTR10i", 2), rep("FOSL_TSSi", 2), rep("GFPi", 2), rep("MCPH1_LTR10i", 2), rep("MEF2D_LTR10i", 2), rep("XRCC4_LTR10i", 2)))
condition

condition = relevel(condition, "GFPi")
condition
# GFPi is now the reference level

# Create a "coldata" table containing the sample names with their appropriate condition
coldata <- data.frame(row.names=colnames(countdata), condition)
coldata

# Construct a DESeqDataSet 
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
dds

# Check that everything looks right
dds$condition

# Double check that GFPi is the control/reference level
dds$condition <- relevel(dds$condition, ref = "GFPi")
dds$condition

# Run deseq2
dds <- DESeq(dds) 
resultsNames(dds)
[1] "Intercept"                      "condition_ATG12_LTR10i_vs_GFPi" "condition_ATG12_TSSi_vs_GFPi"  
[4] "condition_cJun_TSSi_vs_GFPi"    "condition_FGF2_LTR10i_vs_GFPi"  "condition_FOSL_TSSi_vs_GFPi"   
[7] "condition_MCPH1_LTR10i_vs_GFPi" "condition_MEF2D_LTR10i_vs_GFPi" "condition_XRCC4_LTR10i_vs_GFPi"

# Save normalized counts table
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="normalized_counts_crispri.tab", quote = FALSE, row.names = TRUE, sep = "\t")

# Make a pearsons correlation matrix using normalized counts
cormat <- cor(normalized_counts, method="pearson")
cormat
symnum(cormat)

# Plot the correlation matrix
melted_cormat <- melt(cormat)
melted_cormat

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# Make a function to reorder the correlation matrix by clustering
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Now apply the reorder function
cormat <- reorder_cormat(cormat)

# Now replot it, with samples appearing in the clustered order
melted_cormat <- melt(cormat)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# What about a pca plot?
# First, recommended to do a variance stabilizing transformation
# https://master.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#starting-from-count-matrices
vsd <- vst(dds)
class(dds)
class(vsd)
colData(vsd)

# Make a PCA plot (plotPCA is part of the deseq2 package)
nudge <- position_nudge(y = 1.5)
plotPCA(vsd, returnData=FALSE, ntop=20000) + geom_text(aes(label = name), position = nudge) + xlim(-22,28)
