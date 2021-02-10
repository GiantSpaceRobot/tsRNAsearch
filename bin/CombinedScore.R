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
library("xtable")

my.input = read.csv(args[1], sep = "\t", row.names = 1) # Read input

### Convert NAs into 0 (for scores) and 1 (for p-vals)
my.input$slope.score[is.na(my.input$slope.score)] <- 0
my.input$distribution.score[is.na(my.input$distribution.score)] <- 0
my.input$cleavage.score[is.na(my.input$cleavage.score)] <- 0
my.input$Fishers.method.pvalue[is.na(my.input$Fishers.method.pvalue)] <- 1

my.input$Fishers.pval.negLog10 <- -log10(my.input$Fishers.method.pvalue) # Negative log10 of Fisher's method p-values
number.of.genes.under.DESeq2.padj.cutoff <- my.input %>% filter(DESeq2.padj < 0.05) %>% nrow() # Count number of rows with DESeq2 padj under 0.05
### If nothing passed DESeq2 padj threshold, reduce padj values to 0 so DESeq2 does not contribute to combined score
if (number.of.genes.under.DESeq2.padj.cutoff < 1){
  my.input$DESeq2.padj <- 0
  my.input$DESeq2.padj.negLog10 <- 0
} else {
  my.input$DESeq2.padj.negLog10 <- -log10(my.input$DESeq2.padj) # Negative log10 of DESeq2 p-adj
  my.input$DESeq2.padj.negLog10 <- replace(x = my.input$DESeq2.padj.negLog10, is.infinite(my.input$DESeq2.padj.negLog10), 20) # Replace -log10 of (p-adj = 0) with 20 as it is Inf
}
final.df <- my.input %>% select(c(slope.score, distribution.score, cleavage.score, Fishers.pval.negLog10, DESeq2.padj.negLog10)) # Subset DF

### Create intermediate DF to Square root transform Slope, Distribution, and Cleavage scores (14-7-20).
final.df.logs <- final.df
final.df.logs$slope.score <- sqrt(final.df.logs$slope.score) 
final.df.logs$distribution.score <- sqrt(final.df$distribution.score)
final.df.logs$cleavage.score <- sqrt(final.df.logs$cleavage.score)
final.df.logs$DESeq2.padj.negLog10[is.na(final.df.logs$DESeq2.padj.negLog10)] <- 0
#final.df.logs$DESeq2.padj.negLog10 <- final.df.logs$DESeq2.padj.negLog10 + 1 # Add pseudocount

final.percentages.df <- as.data.frame(sapply(final.df.logs, function(x) (x/max(x))*100)) # Calculate percentages for each column
final.percentages.df[is.na(final.percentages.df)] <- 0
rownames(final.percentages.df) <- rownames(final.df.logs) # Take rownames from old DF
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
            file = paste0(args[2], "_relative-score-results.tsv"),
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
pdf(file = paste0(args[2], "_tsRNAs.pdf"), width = 7, height = 5)
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
pdf(file = paste0(args[2], "_ncRNAs.pdf"), width = 7, height = 5)
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
pdf(file = paste0(args[2], "_tRNAs-and-ncRNAs.pdf"), width = 7, height = 5)
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

### Final output DF preparation
original.scores <- my.input
original.scores$feature <- rownames(original.scores) # Add rownames as column
final.percentages.df$feature <- rownames(final.percentages.df) # Add rownames as column
features.and.combined.scores.df <- final.percentages.df %>% select(c(feature, combined.score)) # Subset DF
output.df <- merge(features.and.combined.scores.df, original.scores, by="feature", all.x=TRUE) # Merge original input with combined score DF
rownames(output.df) <- output.df$feature
output.df <- output.df %>% rownames_to_column('genes') %>% 
  select(c(combined.score, slope.score, Fishers.method.pvalue, distribution.score, DESeq2.Log2FC, DESeq2.padj, cleavage.score, genes)) %>%
  arrange(desc(combined.score)) %>%
  column_to_rownames('genes') # Subset new merged DF
output.df$DESeq2.padj[is.na(output.df$DESeq2.padj)] <- 1

### If genes passed DESeq2 padj cutoff and I had to change every value to 0 to generate combined score, change to 1 
if (sum(output.df$DESeq2.padj) == 0){
  output.df$DESeq2.padj = 1
}

# If there are more than 20 features, show top 20
if(nrow(output.df) > 20){
  top.results <- head(output.df, n = 20)
} else {
  top.results <- output.df
}

### Write output
write.table(x = output.df, file = paste0(args[2], "_summarised_absolute-score-results.tsv"), quote = F, sep = "\t", row.names = T)
write.table(x = top.results, file = paste0(args[2], "_summarised_absolute-score-results.tsv"), quote = F, sep = "\t", row.names = T)

### Write top results as HTML (s = string, g = number, f = scientific number)
print(xtable(top.results,              
             display=c("s",   # Feature names 
                       "f",   # Combined score
                       "f",   # Slope score
                       "g",   # Fisher's method p value
                       "f",   # Distribution score
                       "f",   # DESeq2 Log2FC
                       "g",   # DESeq2 padj
                       "f")), # Cleavage score
      math.style.exponents = TRUE, 
      type="html", 
      file = paste0(args[2], "_summarised_top-results.html")) # Write top results to HTML for HTML report generation

