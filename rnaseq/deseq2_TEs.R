# Template R script for DEseq2 analysis on TE-transcripts output (TEs only, no genes)

suppressPackageStartupMessages(library(DESeq2))

# Clear workspace
rm(list = ls())

# Read table
TEdata <- read.table("TNF_24hr_bothReps_TEsonly.counttab",header=T,row.names=1)
head(TEdata)

groups <- factor(c(rep("TGroup",2),rep("CGroup",2)))
min_read <- 1
TEdata <- TEdata[apply(TEdata,1,function(x){max(x)}) > min_read,]

sampleInfo <- data.frame(groups,row.names=colnames(TEdata))

dds <- DESeqDataSetFromMatrix(countData = TEdata, colData = sampleInfo, design = ~groups)

dds$groups = relevel(dds$groups,ref="CGroup")
dds <- DESeq(dds)

# Save the count tables

# Raw counts
raw_counts_tetranscripts_tnf_24hr = counts(dds)
head(raw_counts_tetranscripts_tnf_24hr)
write.table(raw_counts_tetranscripts_tnf_24hr, file="raw_counts_tetranscripts_tnf_24hr.tab", quote = FALSE, row.names = TRUE, sep = "\t")

# Normalized counts
normalized_counts_tetranscripts_tnf_24hr <- counts(dds, normalized=TRUE)
tail(normalized_counts_tetranscripts_tnf_24hr)
write.table(normalized_counts_tetranscripts_tnf_24hr, file="normalized_counts_tetranscripts_tnf_24hr.tab", quote = FALSE, row.names = TRUE, sep = "\t")

# Save the deseq2 output

resultsNames(dds)
#BiocManager::install("apeglm")
#library(apeglm)
deseq2_results_tetranscripts_tnf_vs_ctrl <- lfcShrink(dds, coef="groups_TGroup_vs_CGroup",type="apeglm")
head(deseq2_results_tetranscripts_tnf_vs_ctrl)
dim(deseq2_results_tetranscripts_tnf_vs_ctrl)
#1021

# Remove NAs (e.g. due to no reads in all samples)
deseq2_results_tetranscripts_tnf_vs_ctrl <- na.omit(deseq2_results_tetranscripts_tnf_vs_ctrl)
dim(deseq2_results_tetranscripts_tnf_vs_ctrl)
#1021

# Sort by lowest padj value
deseq2_results_tetranscripts_tnf_vs_ctrl <- deseq2_results_tetranscripts_tnf_vs_ctrl[order(deseq2_results_tetranscripts_tnf_vs_ctrl$padj), ]
head(deseq2_results_tetranscripts_tnf_vs_ctrl, 10)

# Save deseq-normalized table sorted by padj
write.table(deseq2_results_tetranscripts_tnf_vs_ctrl, file="deseq2_results_tetranscripts_tnf_vs_ctrl.tab", quote = FALSE, row.names = TRUE, sep = "\t")
