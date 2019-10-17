#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### Generate text reports
### 
###-------------------------------------------------------------------------------

suppressMessages(library(dplyr))
suppressMessages(library(plyr))

args = commandArgs(trailingOnly=TRUE)

### Load input files
pvals <- read.table(args[1], sep = "\t", header = T)
distribution.scores.tsRNAs <- read.table(args[2], sep = "\t", header = T)
distribution.scores.snomiRNAs <- read.table(args[3], sep = "\t", header = T)
cleavage.scores.tsRNAs <- read.table(args[4], sep = "\t", header = T)
cleavage.scores.snomiRNAs <- read.table(args[5], sep = "\t", header = T)
DEGs <- read.csv(args[6])

colnames(DEGs)[1]  <- "feature" # rename column 1 to features for DEG DF

### Combine tsRNA and snomiRNA dataframes
distribution.scores <- rbind.data.frame(distribution.scores.tsRNAs, distribution.scores.snomiRNAs)
cleavage.scores <- rbind.data.frame(cleavage.scores.tsRNAs, cleavage.scores.snomiRNAs)

### Join dataframes
joined1 <- join(x = pvals, y = distribution.scores, by = "feature")
joined1.subset <- select(joined1, feature, Fishers.method.pvalue, distribution.score)
joined2 <- join(x = cleavage.scores, y = DEGs, by = "feature")
joined2.subset <- select(joined2, feature, "DESeq2.Log2FC" = log2FoldChange, "DESeq2.pvalue" = pvalue, "DESeq2.padj" = padj, cleavage.score)
final.df <- join(x = joined1.subset, y = joined2.subset, by = "feature")
final.df$distribution.score <- as.numeric(final.df$distribution.score) # Convert distribution.score column to numeric
final.df <- final.df[order(-final.df$distribution.score),] # Sort by distribution score

# If there are more than 20 features, show top 20
if(nrow(final.df) > 20){
  top.results <- head(final.df, n = 20)
} else {
  top.results <- final.df
}

### Playing around with data relationships
#plot(x = abs(final.df$DESeq2.Log2FC), -log10(final.df$DESeq2.pvalue))
#pairs(final.df)

write.table(x = final.df, file = paste0(args[7], ".summarised.all-results.txt"), quote = F, sep = "\t", row.names = F)
write.table(x = top.results, file = paste0(args[7], ".summarised.top-results.txt"), quote = F, sep = "\t", row.names = F)
