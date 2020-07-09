#!/usr/bin/env Rscript

## Script name:CombinedScore.R
##
## Purpose of script: Generate combined score
##
## Author: Dr Paul Donovan
##
## Date Created: 9-7-20
##
## Copyright: MIT license
##
## Email: pauldonovan@rcsi.com

args = commandArgs(trailingOnly=TRUE)

# args 1 = Score input file
# args 2 = Output file prefix

library("dplyr")
library("ggplot2")
library("tibble")

my.input = read.csv(args[1], sep = "\t", row.names = 1) # Read input
#my.input <- read.csv("/home/paul/Documents/Pipelines/Analyses_tsRNAsearch/Test_26-6-20/Data/CytC_vs_TotalRNA.summarised.all-results.txt", 
#                     sep = "\t", row.names = 1)
my.input[is.na(my.input)] <- 0 # Convert NAs to 0
my.input$Fishers.pval.negLog10 <- -log10(my.input$Fishers.method.pvalue) # Negative log10 of Fisher's method p-values
my.input$DESeq2.padj.negLog10 <- -log10(my.input$DESeq2.padj) # Negative log10 of DESeq2 p-adj
my.input$DESeq2.padj.negLog10 <- replace(x = my.input$DESeq2.padj.negLog10, is.infinite(my.input$DESeq2.padj.negLog10), 20) # Replace -log10 of (p-adj = 0) with 20 as it is Inf
final.df <- my.input %>% select(c(slope.score, distribution.score, cleavage.score, Fishers.pval.negLog10, DESeq2.padj.negLog10)) # Subset DF

final.percentages.df <- as.data.frame(sapply(final.df, function(x) (x/max(x))*100)) # Calculate percentages for each column
rownames(final.percentages.df) <- rownames(final.df) # Take rownames from old DF
final.percentages.df$combined.score <- rowSums(final.percentages.df) # Create combined scores
final.percentages.df <- final.percentages.df[,c(6,1:5)]
final.percentages.df <- final.percentages.df %>% 
  rownames_to_column('genes') %>% 
  arrange(desc(combined.score)) %>%
  column_to_rownames('genes') # Sort in descending order
is.num <- sapply(final.percentages.df, is.numeric)
final.percentages.df[is.num] <- lapply(final.percentages.df[is.num], round, 2) # Round all numbers to 2 decimal places

# Write output
write.table(x = final.percentages.df, 
            file = paste0(args[2], ".all-results.tsv"),
            quote = FALSE, 
            sep = "\t", 
            col.names=NA)

### Prepare for plot

# If there are more than 20 features, show top 20 (of tsRNAs)
tsRNA.subset <- final.percentages.df[!grepl("ENS", rownames(final.percentages.df)),]
if(nrow(tsRNA.subset) > 20){
  newdata.subset <- head(tsRNA.subset, n = 20)
} else {
  newdata.subset <- tsRNA.subset
}
newdata.subset$feature <- rownames(newdata.subset)
newdata.subset$feature <- factor(newdata.subset$feature, levels = newdata.subset$feature[order(newdata.subset$combined.score)]) # Refactorise based on order
pdf(file = paste0(args[2], ".tsRNAs.pdf"), width = 7, height = 5)
ggplot(data = newdata.subset, mapping = aes(feature, 
                                            newdata.subset$combined.score, 
                                            color=newdata.subset$combined.score)) +
  geom_point() +
  ggtitle("Combined score") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size=10)) +
  scale_color_gradient(low="blue", high="red") +
  scale_x_discrete(limits = (levels(newdata.subset$combined.score))) +
  coord_flip() +
  labs(colour = "Combined\n    score", 
       x = "tRNA", 
       y = "Combined score", 
       subtitle = "Max number of features shown = 20")
dev.off()

# If there are more than 20 features, show top 20 (of tsRNAs)
ncRNA.subset <- final.percentages.df[grep("ENS", rownames(final.percentages.df)),]
if(nrow(ncRNA.subset) > 20){
  newdata.subset <- head(ncRNA.subset, n = 20)
} else {
  newdata.subset <- ncRNA.subset
}
newdata.subset$feature <- rownames(newdata.subset)
newdata.subset$feature <- factor(newdata.subset$feature, levels = newdata.subset$feature[order(newdata.subset$combined.score)]) # Refactorise based on order
pdf(file = paste0(args[2], ".ncRNAs.pdf"), width = 7, height = 5)
ggplot(data = newdata.subset, mapping = aes(feature, 
                                            newdata.subset$combined.score, 
                                            color=newdata.subset$combined.score)) +
  geom_point() +
  ggtitle("Combined score") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size=10)) +
  scale_color_gradient(low="blue", high="red") +
  scale_x_discrete(limits = (levels(newdata.subset$combined.score))) +
  coord_flip() +
  labs(colour = "Combined\n    score", 
       x = "ncRNA", 
       y = "Combined score", 
       subtitle = "Max number of features shown = 20")
dev.off()

# If there are more than 20 features, show top 20
if(nrow(final.percentages.df) > 20){
  newdata.subset <- head(final.percentages.df, n = 20)
} else {
  newdata.subset <- final.percentages.df
}
newdata.subset$feature <- rownames(newdata.subset)
newdata.subset$feature <- factor(newdata.subset$feature, levels = newdata.subset$feature[order(newdata.subset$combined.score)]) # Refactorise based on order
pdf(file = paste0(args[2], ".tRNAs-and-ncRNAs.pdf"), width = 7, height = 5)
ggplot(data = newdata.subset, mapping = aes(feature, 
                                            newdata.subset$combined.score, 
                                            color=newdata.subset$combined.score)) +
  geom_point() +
  ggtitle("Combined score") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size=8)) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0, size=8)) +
  scale_color_gradient(low="blue", high="red") +
  scale_x_discrete(limits = (levels(newdata.subset$combined.score))) +
  coord_flip() +
  labs(colour = "Combined\n    score", 
       x = "ncRNA", 
       y = "Combined score", 
       subtitle = "Max number of features shown = 20")
dev.off()
