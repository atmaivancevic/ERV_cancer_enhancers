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
write.table(raw_counts_24hr, file="/scratch/Users/ativ2716/data/1_testing_github_code/raw_counts_24hr.tab", quote = FALSE, row.names = TRUE, sep = "\t")

# Normalized counts
normalized_counts_24hr <- counts(dds, normalized=TRUE)
tail(normalized_counts_24hr)
write.table(normalized_counts_24hr, file="/scratch/Users/ativ2716/data/1_testing_github_code/normalized_counts_24hr.tab", quote = FALSE, row.names = TRUE, sep = "\t")
