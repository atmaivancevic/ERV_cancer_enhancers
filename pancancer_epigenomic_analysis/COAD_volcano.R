### Install and load DESeq2.
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#library("DESeq2")
#biocLite("pheatmap")
#library("pheatmap")
#library("ggplot2")
#biocLite("ggrepel")
#library("ggrepel")

### Read in the output from featureCounts.
countdata <- read.table("/scratch/Users/ativ2716/data/TCGA_subtract_Roadmap/cancer_minus_roadmap167categories/COAD_minus_Roadmap167.tab.volcano", header=TRUE, sep="\t")
head(countdata)

# Set first column to be the rownames
rownames(countdata) <- countdata[,1]
head(countdata)

### Remove the first column
countdata <- countdata[ ,2:ncol(countdata)]
head(countdata)

# volcano plot
fishers_two_tail_pval <- countdata[,c(2)]
head(fishers_two_tail_pval)
hist(fishers_two_tail_pval)
hist(-log10(fishers_two_tail_pval))

odds_ratio <- countdata[,c(1)]
head(odds_ratio)
hist(odds_ratio)

new <- countdata[,c(1,2)]
head(new)
new = as.data.frame(new)
new$probename <- rownames(new)
head(new)

new["LTR10F", ]
new["LTR10A", ]

# Separate data into subsets
value1 = subset(new, fishers_two_tail_pval<0.05)
nrow(value1)
head(value1)

value2 = subset(new, odds_ratio>1)
nrow(value2)
value2

value3 = subset(new, fishers_two_tail_pval<0.05 & odds_ratio>1)
dim(value3)
value3

# this last subset will be for labelling the most significant genes
value4 = subset(new, fishers_two_tail_pval<0.00000000000000000000000005 & odds_ratio>10)
nrow(value4)
rownames(value4)

value5 = subset(new, fishers_two_tail_pval<0.0000000000000000000005 & odds_ratio<0.25)
nrow(value5)
rownames(value5)

# plot it
ggplot(new, aes(log2(odds_ratio), -log10(fishers_two_tail_pval)), colour="grey") +
  scale_color_discrete(name = 'Labels') +
  theme_classic() +
  labs(y="-Log10 P-value", x = "Log2 Odds Ratio") +
  # Set all dots color to grey
  geom_point(data=new, colour = "grey", size = 1) +
  # If pvalue<0.05, change dot color to green
  geom_point(data=new[which(new $odds_ratio >3 & new $fishers_two_tail_pval <0.005),], colour = "red", size=1) +
  geom_point(data=new[which(new $odds_ratio <(1/3) & new $fishers_two_tail_pval <0.005),], colour = "blue", size=1) +
  geom_text_repel(data = value4, mapping = aes(log2(odds_ratio), -log10(fishers_two_tail_pval), label = rownames(value4)), force = 1) +
  geom_text_repel(data = value5, mapping = aes(log2(odds_ratio), -log10(fishers_two_tail_pval), label = rownames(value5)),force = 1) +
  xlim(-8.2, 8.2) +
  #ylim(0, 66) +
  theme(text = element_text(size=11, family="Helvetica"), 
        axis.text.x = element_text(size = 10, family="Helvetica"),
        axis.text.y = element_text(size = 10, family="Helvetica"),
        axis.title.x = element_text(size = 10, family="Helvetica"),
        axis.title.y = element_text(size = 10, family="Helvetica"))
