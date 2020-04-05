#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### Generate text reports
### 
###-------------------------------------------------------------------------------

suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(xtable))

args = commandArgs(trailingOnly=TRUE)

### Load input files
pvals <- read.table(args[1], sep = "\t", header = T)
distribution.scores.tsRNAs <- read.table(args[2], sep = "\t", header = T)
distribution.scores.ncRNAs <- read.table(args[3], sep = "\t", header = T)
cleavage.scores.tsRNAs <- read.table(args[4], sep = "\t", header = T)
cleavage.scores.ncRNAs <- read.table(args[5], sep = "\t", header = T)
DEGs <- read.csv(args[6])
slope.scores.tsRNAs <- read.table(args[7], sep = "\t", header = T)
slope.scores.ncRNAs <- read.table(args[8], sep = "\t", header = T)

# pvals <- read.table("/home/paul/Documents/Pipelines/Analyses_tsRNAsearch/PathTest/Data/FishersMethod_results/Feature-P-values_FisherMethod_pvalues.tsv", sep = "\t", header = T)
# distribution.scores.tsRNAs <- read.table("/home/paul/Documents/Pipelines/Analyses_tsRNAsearch/PathTest/Data/Distribution_results/CytC_vs_TotalRNA_High-distribution-tsRNAsTEST.txt", sep = "\t", header = T)
# distribution.scores.ncRNAs <- read.table("/home/paul/Documents/Pipelines/Analyses_tsRNAsearch/PathTest/Data/Distribution_results/CytC_vs_TotalRNA_High-distribution-ncRNAsTEST.txt", sep = "\t", header = T)
# cleavage.scores.tsRNAs <- read.table("/home/paul/Documents/Pipelines/Analyses_tsRNAsearch/PathTest/Data/Cleavage_results/CytC_vs_TotalRNA_High-Cleavage-tRNAs.tsv", sep = "\t", header = T)
# cleavage.scores.ncRNAs <- read.table("/home/paul/Documents/Pipelines/Analyses_tsRNAsearch/PathTest/Data/Cleavage_results/CytC_vs_TotalRNA_High-Cleavage-ncRNAs.tsv", sep = "\t", header = T)
# DEGs <- read.csv("/home/paul/Documents/Pipelines/Analyses_tsRNAsearch/PathTest/Data/DE_Results/DESeq2/CytC_vs_TotalRNA_DESeq2-output.csv")
# slope.scores.tsRNAs <- read.table("/home/paul/Documents/Pipelines/Analyses_tsRNAsearch/PathTest/Data/Slope_results/CytC_vs_TotalRNA_High-Slope-tRNAs.tsv", sep = "\t", header = T)
# slope.scores.ncRNAs <- read.table("/home/paul/Documents/Pipelines/Analyses_tsRNAsearch/PathTest/Data/Slope_results/CytC_vs_TotalRNA_High-Slope-ncRNAs.tsv", sep = "\t", header = T)

colnames(DEGs)[1]  <- "feature" # rename column 1 to features for DEG DF

### Combine tsRNA and ncRNA dataframes
distribution.scores <- rbind.data.frame(distribution.scores.tsRNAs, distribution.scores.ncRNAs)
cleavage.scores <- rbind.data.frame(cleavage.scores.tsRNAs, cleavage.scores.ncRNAs)
slope.scores <- rbind.data.frame(slope.scores.tsRNAs, slope.scores.ncRNAs)
#DEGs.cutoff <- DEGs %>% filter(padj < 0.05)

### Join dataframes
joined1 <- join(x = pvals, y = distribution.scores, by = "feature")
joined1.subset <- select(joined1, feature, Fishers.method.pvalue, distribution.score)
joined2 <- join(x = cleavage.scores, y = DEGs, by = "feature")
joined2.subset <- select(joined2, feature, "DESeq2.Log2FC" = log2FoldChange, "DESeq2.pvalue" = pvalue, "DESeq2.padj" = padj, cleavage.score)
joined3 <- join(x = slope.scores, y = joined1.subset, by = "feature")
joined3.subset <- select(joined3, feature, Fishers.method.pvalue, distribution.score, slope.score)
final.df <- join(x = joined3.subset, y = joined2.subset, by = "feature")
final.df$distribution.score <- as.numeric(final.df$distribution.score) # Convert distribution.score column to numeric
# Add rank columns for each of the metrics:
#final.df <- final.df %>% mutate(rank.distribution = dense_rank(desc(distribution.score)),
#                                rank.cleavage = dense_rank(desc(cleavage.score)),
#                                rank.slope = dense_rank(desc(slope.score)),
#                                rank.fisher = dense_rank(desc(Fishers.method.pvalue)))#,
                                #rank.DE = dense_rank(desc(DESeq2.padj)))

### Fill NAs with poor ranks so features with a score in only one method get good metascore
### First, try replace NAs with no. of features found using each method
### If no features were found using this method, fill with total no. of features found
# Get numbers from each method
# distribution.number.of.features <- ifelse((nrow(distribution.scores)) > 0, 
#                                           nrow(distribution.scores),
#                                           nrow(final.df))
# slope.number.of.features <- ifelse((nrow(slope.scores)) > 0, 
#                                           nrow(slope.scores),
#                                           nrow(final.df))
# cleavage.number.of.features <- ifelse((nrow(cleavage.scores)) > 0, 
#                                           nrow(cleavage.scores),
#                                           nrow(final.df))
#DE.number.of.features <- ifelse((nrow(DEGs.cutoff)) > 0, 
#                                          nrow(DEGs.cutoff),
#                                          nrow(final.df))
# pvals.number.of.features <- ifelse((nrow(pvals)) > 0, 
#                                           nrow(pvals),
#                                           nrow(final.df))
# Replacement step
# final.df$rank.distribution[is.na(final.df$rank.distribution)] <- distribution.number.of.features
# final.df$rank.cleavage[is.na(final.df$rank.cleavage)] <- cleavage.number.of.features
# final.df$rank.slope[is.na(final.df$rank.slope)] <- slope.number.of.features
# final.df$rank.fisher[is.na(final.df$rank.fisher)] <- pvals.number.of.features
#final.df$rank.DE[is.na(final.df$rank.DE)] <- DE.number.of.features

# Get the mean of ranks and create new column called metascore:
#final.df <- final.df %>% mutate(metascore = rowMeans(select(final.df, starts_with("rank")), na.rm = TRUE))
#final.df <- final.df %>% select(feature, metascore, everything()) #Bring metascore to start of dataframe
#final.df <- final.df[order(final.df$metascore),] # Sort by metascore
final.df <- final.df %>% select(feature, slope.score, everything()) #Bring slope.score to start of dataframe
final.df <- final.df[order(-final.df$slope.score),] # Sort by metascore
write.table(x = final.df, file = paste0(args[9], ".summarised.all-results.all-columns.tsv"), quote = F, sep = "\t", row.names = F)
#final.df <- final.df %>% select(-starts_with("rank")) # Remove all rank columns
#final.df <- final.df[order(-final.df$distribution.score),] # Sort by distribution score

# If there are more than 20 features, show top 20
if(nrow(final.df) > 20){
  top.results <- head(final.df, n = 20)
} else {
  top.results <- final.df
}

write.table(x = final.df, file = paste0(args[9], ".summarised.all-results.txt"), quote = F, sep = "\t", row.names = F)
write.table(x = top.results, file = paste0(args[9], ".summarised.top-results.txt"), quote = F, sep = "\t", row.names = F)
rownames(top.results) <- top.results$feature # Replace rownames with feature names
top.results <- top.results %>% select(-feature) # Remove feature names (column 1)

### Write top results as HTML (s = string, g = number, f = scientific number)
print(xtable(top.results,              
             display=c("s",   # Feature names 
                       "f",   # Slope score
                       "g",   # Fisher's method p value
                       "f",   # Distribution score
                       "f",   # DESeq2 Log2FC
                       "g",   # DESeq2 p-value
                       "g",   # DESeq2 padj
                       "f")), # Cleavage score
      math.style.exponents = TRUE, 
      type="html", 
      file = paste0(args[9], ".summarised.top-results.html")) # Write top results to HTML for HTML report generation
