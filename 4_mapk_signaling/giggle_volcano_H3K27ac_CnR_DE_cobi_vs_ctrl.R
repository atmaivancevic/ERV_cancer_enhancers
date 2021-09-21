library("ggplot2")
library("ggrepel")

# Read in the table
giggle_CnR_cobi <- read.table("DE_regions_H3K27ac_CnR_cobimetinib_giggle_results.tab", header=TRUE, sep="\t", row.names=1)
head(giggle_CnR_cobi)

# Fix the column names that got right-shifted
colnames(giggle_CnR_cobi) = c("repeat_family_size", "overlaps", "odds_ratio", "fishers_two_tail", "fishers_left_tail", "fishers_right_tail", "giggle_score")
giggle_CnR_cobi[,ncol(giggle_CnR_cobi)] <- NULL
head(giggle_CnR_cobi)

# Subset the most enriched/depleted TEs (for labelling on the plot)
labelsEnriched = subset(giggle_CnR_cobi, fishers_two_tail<0.00000005 & odds_ratio>9)
nrow(labelsEnriched)
rownames(labelsEnriched)

labelsDepleted = subset(giggle_CnR_cobi, fishers_two_tail<0.0000005 & odds_ratio<0.1)
nrow(labelsDepleted)
rownames(labelsDepleted)

# Plot volcano
ggplot(giggle_CnR_cobi, aes(log2(odds_ratio),-log10(fishers_two_tail)), colour="grey") +
  scale_color_discrete(name='Labels') +
  theme_classic() +
  labs(y="-Log10 P-value", x="Log2 Odds Ratio") +
  ggtitle("TEs enriched in H3K27ac-marked C&R DE regions (Cobimetinib)") +
  # Set all dots to grey
  geom_point(data=giggle_CnR_cobi, colour="grey", size=1) +
  # Change dot colour to red (enriched) or blue (depleted) based on odds ratio and pval cutoffs
  geom_point(data=giggle_CnR_cobi[which(giggle_CnR_cobi$odds_ratio>3 & giggle_CnR_cobi$fishers_two_tail<0.05),], colour="red", size=1) +
  geom_point(data=giggle_CnR_cobi[which(giggle_CnR_cobi$odds_ratio<(1/3) & giggle_CnR_cobi$fishers_two_tail<0.05),], colour="blue", size=1) +
  # Add labels to the most enriched/depleted points
  geom_text_repel(data=labelsEnriched, mapping=aes(log2(odds_ratio),-log10(fishers_two_tail),label=rownames(labelsEnriched)), force=1) +
  geom_text_repel(data=labelsDepleted, mapping=aes(log2(odds_ratio),-log10(fishers_two_tail),label=rownames(labelsDepleted)), force=1) +
  xlim(-9.2,9.2) +
  theme(text=element_text(size=11,family="Helvetica"), 
        axis.text.x=element_text(size=10,family="Helvetica"),
        axis.text.y=element_text(size=10,family="Helvetica"),
        axis.title.x=element_text(size=10,family="Helvetica"),
        axis.title.y=element_text(size=10,family="Helvetica"))
