library("ggplot2")
library("ggrepel")

# Read in the table
COAD_minus_Roadmap <- read.table("COAD_minus_Roadmap.tab.volcano", header=TRUE, sep="\t")
head(COAD_minus_Roadmap)

# Set first column to be the rownames
rownames(COAD_minus_Roadmap) <- COAD_minus_Roadmap[,1]
head(COAD_minus_Roadmap)

# Remove the first column
COAD_minus_Roadmap <- COAD_minus_Roadmap[ ,2:ncol(COAD_minus_Roadmap)]
head(COAD_minus_Roadmap)

# Plot histograms & volcano plot
fishers_two_tail_pval <- COAD_minus_Roadmap[,c(2)]
head(fishers_two_tail_pval)
hist(fishers_two_tail_pval)
hist(-log10(fishers_two_tail_pval))

odds_ratio <- COAD_minus_Roadmap[,c(1)]
hist(odds_ratio)

COAD_simple <- COAD_minus_Roadmap[,c(1,2)]
head(COAD_simple)
COAD_simple = as.data.frame(COAD_simple)
COAD_simple$probename <- rownames(COAD_simple)
head(COAD_simple)

# This last subset will be for adding labels to the most enriched/depleted TEs
labelsEnriched = subset(COAD_simple, fishers_two_tail_pval<0.00000000000000000000000005 & odds_ratio>30)
nrow(labelsEnriched)
rownames(labelsEnriched)

labelsDepleted = subset(COAD_simple, fishers_two_tail_pval<0.0000000000000000000005 & odds_ratio<0.25)
nrow(labelsDepleted)
rownames(labelsDepleted)

# Plot it
ggplot(COAD_simple, aes(log2(odds_ratio),-log10(fishers_two_tail_pval)), colour="grey") +
  scale_color_discrete(name='Labels') +
  theme_classic() +
  labs(y="-Log10 P-value", x="Log2 Odds Ratio") +
  # Set all dots to grey
  geom_point(data=COAD_simple, colour="grey", size=1) +
  # Change dot colour to red (enriched) or blue (depleted) based on odds ratio and pval cutoffs
  geom_point(data=COAD_simple[which(COAD_simple$odds_ratio>3 & COAD_simple$fishers_two_tail_pval<0.005),], colour="red", size=1) +
  geom_point(data=COAD_simple[which(COAD_simple$odds_ratio<(1/3) & COAD_simple$fishers_two_tail_pval<0.005),], colour="blue", size=1) +
  # Add labels to the most enriched/depleted points
  geom_text_repel(data=labelsEnriched, mapping=aes(log2(odds_ratio),-log10(fishers_two_tail_pval),label=rownames(labelsEnriched)), force=1) +
  geom_text_repel(data=labelsDepleted, mapping=aes(log2(odds_ratio),-log10(fishers_two_tail_pval),label=rownames(labelsDepleted)), force=1) +
  xlim(-8.2,8.2) +
  theme(text=element_text(size=11,family="Helvetica"), 
        axis.text.x=element_text(size=10,family="Helvetica"),
        axis.text.y=element_text(size=10,family="Helvetica"),
        axis.title.x=element_text(size=10,family="Helvetica"),
        axis.title.y=element_text(size=10,family="Helvetica"))
