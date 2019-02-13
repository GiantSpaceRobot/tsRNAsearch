#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script calculates statistics for 5' vs 3' read coverage distributions
### 
###-------------------------------------------------------------------------------

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
} 

#### Input file
input <- read.table(args[1])
#input <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/subset/Results/Ang_1_ATCACG_L008_R1_001/snomiRNA-alignment/accepted_hits_sorted.depth")

df <- split( input , f = input$V1 )  # Split dataframe based on column 1 elements

results.df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("feature",
                                                                 "Mean of 5prime",
                                                                 "Mean of 3prime",
                                                                 #"Wilcox p-value",
                                                                 "Kolmogorov-smirnov p-value",
                                                                 "-Log10 of KS p-value")) # Initialise empty dataframe with column headers

for(subset in df) {
  feature <- as.character(subset[1,1])
  subset.length <- nrow(subset)
  half.length <- as.integer(subset.length/2)
  if (half.length < 150) {
    length.penalty <- half.length/150###
  } else {
    length.penalty <- 1
  }
  fiveprime <- subset[1:half.length,]
  five.distribution <- fiveprime$V3
  threeprime <- subset[(half.length+1):subset.length,]
  three.distribution <- threeprime$V3
  mean1 <- mean(threeprime$V3)
  mean2 <- mean(fiveprime$V3)
  ### If neither mean1 nor mean2 have a mean coverage over 10, default the KS p-value to 0
  if (mean1 > 10 | mean2 > 10){
    ks.output <- ks.test(five.distribution, three.distribution)
    ks.pvalue <- ks.output$p.value
    #wilcox.output <- wilcox.test(five.distribution, three.distribution)
    #wilcox.pvalue <- wilcox.output$p.value
    if (ks.pvalue == 0) {   # If the ks p-value is 0 (and subsequently the -Log10 is also 0), artificially insert a very low value
      ks.pvalue <- 1e-20
    }
    #minusLogKS <- -log10(ks.pvalue)
    minusLogKS <- (-log10(ks.pvalue))*length.penalty
  } else {
    wilcox.pvalue <- "NA"
    ks.pvalue <- 0
    minusLogKS <- 0
  }
  results.df[nrow(results.df) + 1,] = list(feature, mean1, mean2, ks.pvalue, minusLogKS)

}

results.df <- results.df[order(-results.df$`-Log10 of KS p-value`),]
#newdata <- results.df[ which(results.df$`-Log10 of KS p-value` > 1.30103), ]
newdata <- results.df[ which(results.df$`-Log10 of KS p-value` > 3), ] # Equivalent to p-value of 0.001. More stringent than 0.05


write.table(results.df, 
            file = paste0(args[2], ".all-features.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
write.table(newdata, 
            file = paste0(args[2], ".potentially-cleaved-features.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

### Make plots
#png.width <- nrow(results.df)*15
#png(file = paste0(args[2], ".all-features.png"), width = png.width, height = 400)
#ggplot(data = results.df, mapping = aes(feature, `-Log10 of KS p-value`, color=`-Log10 of KS p-value`)) +
#  geom_point() +
#  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#  scale_color_gradient(low="blue", high="red")
#dev.off()

pdf.width <- nrow(newdata)*0.2 + 3
pdf(file = paste0(args[2], ".potentially-cleaved-features.pdf"), width = pdf.width, height = 5)
ggplot(data = newdata, mapping = aes(feature, `-Log10 of KS p-value`, color=`-Log10 of KS p-value`)) +
  geom_point() +
  ggtitle("Feature cleavage likelihood (-Log10 cutoff = 3)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_gradient(low="blue", high="red")
dev.off()


