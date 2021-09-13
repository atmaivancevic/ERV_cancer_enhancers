# Clear workspace
rm(list = ls())

# read table
countdata <- read.table("featureCounts_withGeneNames.txt", sep="\t", header = TRUE)
head(countdata)

# Set second column to be the rownames
rownames(countdata) <- countdata[,2]
head(countdata)

# Remove the first four columns (geneid, chr, description, length)
countdata <- countdata[ ,5:ncol(countdata)]
head(countdata)

# Convert countdata table into a matrix - necessary for running DESeq2.
countdata <- as.matrix(countdata)
head(countdata)

# Assign control vs test samples
condition <- factor(c(rep("cobi_24", 2), rep("tnf_24", 2), rep("ctrl_24", 2)))
condition

# Relevel to make ctrl the reference level
condition = relevel(condition, "ctrl_24")
condition

# Create a "coldata" table containing the sample names with their appropriate condition
coldata <- data.frame(row.names=colnames(countdata), condition)
coldata

# Construct a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
dds

# Check that everything looks right
dds$condition

dds <- DESeq(dds)
resultsNames(dds)

# Save the count tables

# Raw counts
raw_counts_24hr = counts(dds)
head(raw_counts_24hr)
write.table(raw_counts_24hr, file="raw_counts_24hr.tab", quote = FALSE, row.names = TRUE, sep = "\t")

# Normalized counts
normalized_counts_24hr <- counts(dds, normalized=TRUE)
tail(normalized_counts_24hr)
write.table(normalized_counts_24hr, file="normalized_counts_24hr.tab", quote = FALSE, row.names = TRUE, sep = "\t")

###########################################################################################################

# 1. Compare Cobimetinib to control (untreated)

resultsNames(dds)

#BiocManager::install("apeglm")
#library(apeglm)
deseq2_results_cobi_vs_ctrl <- lfcShrink(dds, coef="condition_cobi_24_vs_ctrl_24",type="apeglm")
head(deseq2_results_cobi_vs_ctrl)

# Remove NAs (e.g. due to no reads in all samples)
deseq2_results_cobi_vs_ctrl <- na.omit(deseq2_results_cobi_vs_ctrl)
head(deseq2_results_cobi_vs_ctrl)
dim(deseq2_results_cobi_vs_ctrl)
#16553

# Sort by lowest padj value
deseq2_results_cobi_vs_ctrl <- deseq2_results_cobi_vs_ctrl[order(deseq2_results_cobi_vs_ctrl$padj),]
head(deseq2_results_cobi_vs_ctrl, 20)

# Save deseq-normalized table sorted by padj
write.table(deseq2_results_cobi_vs_ctrl, file="deseq2_results_cobi_vs_ctrl.tab", quote = FALSE, row.names = TRUE, sep = "\t")

###########################################################################################################

# 2. Compare TNF-alpha to control (untreated)

resultsNames(dds)
deseq2_results_tnf_vs_ctrl <- lfcShrink(dds, coef="condition_tnf_24_vs_ctrl_24",type="apeglm")
head(deseq2_results_tnf_vs_ctrl)

# Remove NAs (e.g. due to no reads in all samples)
deseq2_results_tnf_vs_ctrl <- na.omit(deseq2_results_tnf_vs_ctrl)
head(deseq2_results_tnf_vs_ctrl)
dim(deseq2_results_tnf_vs_ctrl)
#14790

# Sort by lowest padj value
deseq2_results_tnf_vs_ctrl <- deseq2_results_tnf_vs_ctrl[order(deseq2_results_tnf_vs_ctrl$padj), ]
head(deseq2_results_tnf_vs_ctrl, 10)

# Save deseq-normalized table sorted by padj
write.table(deseq2_results_tnf_vs_ctrl, file="deseq2_results_tnf_vs_ctrl.tab", quote = FALSE, row.names = TRUE, sep = "\t")

# Visualize comparisons using MA plots, count plots, volcanos
# See e.g. Figure 3 R-scripts
