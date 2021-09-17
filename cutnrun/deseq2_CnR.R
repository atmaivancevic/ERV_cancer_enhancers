# Tempate script for running DEseq2 on H3K27ac Cut&Run regions

library("DESeq2")
library("ggplot2")
library("ggrepel")
library("apeglm")

# Read in the tab-delimited count table
K27ac_CnR <- read.table("H3K27ac_CnR_bam_counts_within_peaks.tab", sep="\t", header = FALSE)

# Check the table looks correct
head(K27ac_CnR)

# Create a new column which is the chr and coordinates combined
K27ac_CnR$chr <- paste(K27ac_CnR$V1,K27ac_CnR$V2, sep=":")
K27ac_CnR$coord <- paste(K27ac_CnR$chr,K27ac_CnR$V3, sep="-")

# Set this to be the rownames
rownames(K27ac_CnR) <- K27ac_CnR$coord
head(K27ac_CnR)

# Remove the first four columns (chr, start, stop, regionlabel)
# And the last two columns (coord, chr) since we don't need them anymore
# Retaining only the counts
K27ac_CnR <- K27ac_CnR[, c(5, 6, 7, 8, 9, 10)] 
head(K27ac_CnR)

# Add bam names as column names
colnames(K27ac_CnR) <- c("UT_K27ac_rep1","UT_K27ac_rep2","Cobi_K27ac_rep1","Cobi_K27ac_rep2","TNFa_K27ac_rep1","TNFa_K27ac_rep2")
head(K27ac_CnR)

# Convert table to matrix format
K27ac_CnR <- as.matrix(K27ac_CnR)
head(K27ac_CnR)

# Assign control vs treated samples
condition <- factor(c(rep("untreated", 2), rep("cobimetinib", 2), rep("tnfalpha", 2)))

# Create a "coldata" table containing the sample names with their condition
coldata <- data.frame(row.names=colnames(K27ac_CnR), condition)
coldata

# Construct a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = K27ac_CnR, colData = coldata, design = ~condition)
dds

# Check that everything looks right
dds$condition

# Relevel to set the controls as the reference levels
dds$condition <- relevel(dds$condition, ref = "untreated")
dds$condition

# Run deseq2
dds <- DESeq(dds)
resultsNames(dds)

# Save some tables
# Raw counts 
raw_counts_H3K27ac_CnR = counts(dds)
head(raw_counts_H3K27ac_CnR)
write.table(raw_counts_H3K27ac_CnR, file="raw_counts_H3K27ac_CnR.tab", quote = FALSE, row.names = TRUE, sep = "\t")

# Normalized counts
normalized_counts_H3K27ac_CnR <- counts(dds, normalized=TRUE)
write.table(normalized_counts_H3K27ac_CnR, file="normalized_counts_H3K27ac_CnR.tab", quote = FALSE, row.names = TRUE, sep = "\t")

###########################################################################################################

# 1. Compare Cobimetinib to control (untreated)

resultsNames(dds)

deseq2_results_H3K27ac_CnR_cobi_vs_ctrl <- lfcShrink(dds, coef="condition_cobimetinib_vs_untreated",type="apeglm")
head(deseq2_results_H3K27ac_CnR_cobi_vs_ctrl)

# Remove NAs (e.g. due to no reads in all samples)
deseq2_results_H3K27ac_CnR_cobi_vs_ctrl <- na.omit(deseq2_results_H3K27ac_CnR_cobi_vs_ctrl)
head(deseq2_results_H3K27ac_CnR_cobi_vs_ctrl)
dim(deseq2_results_H3K27ac_CnR_cobi_vs_ctrl)
#29297

# Sort by lowest padj value
deseq2_results_H3K27ac_CnR_cobi_vs_ctrl <- deseq2_results_H3K27ac_CnR_cobi_vs_ctrl[order(deseq2_results_H3K27ac_CnR_cobi_vs_ctrl$padj),]
head(deseq2_results_H3K27ac_CnR_cobi_vs_ctrl, 20)

# Save deseq-normalized table sorted by padj
write.table(deseq2_results_H3K27ac_CnR_cobi_vs_ctrl, file="deseq2_results_H3K27ac_CnR_cobi_vs_ctrl.tab", quote = FALSE, row.names = TRUE, sep = "\t")

###########################################################################################################

# 2. Compare TNF-alpha to control (untreated)

resultsNames(dds)

deseq2_results_H3K27ac_CnR_tnf_vs_ctrl <- lfcShrink(dds, coef="condition_tnfalpha_vs_untreated",type="apeglm")
head(deseq2_results_H3K27ac_CnR_tnf_vs_ctrl)

# Remove NAs (e.g. due to no reads in all samples)
deseq2_results_H3K27ac_CnR_tnf_vs_ctrl <- na.omit(deseq2_results_H3K27ac_CnR_tnf_vs_ctrl)
head(deseq2_results_H3K27ac_CnR_tnf_vs_ctrl)
dim(deseq2_results_H3K27ac_CnR_tnf_vs_ctrl)
#26865

# Sort by lowest padj value
deseq2_results_H3K27ac_CnR_tnf_vs_ctrl <- deseq2_results_H3K27ac_CnR_tnf_vs_ctrl[order(deseq2_results_H3K27ac_CnR_tnf_vs_ctrl$padj),]
head(deseq2_results_H3K27ac_CnR_tnf_vs_ctrl, 20)

# Save deseq-normalized table sorted by padj
write.table(deseq2_results_H3K27ac_CnR_tnf_vs_ctrl, file="deseq2_results_H3K27ac_CnR_tnf_vs_ctrl.tab", quote = FALSE, row.names = TRUE, sep = "\t")
