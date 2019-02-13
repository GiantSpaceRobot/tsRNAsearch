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
input1 <- read.table(args[1])
input2 <- read.table(args[2])

#input1 <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/subset/Results/Data/Intermediate-files/tiRNA.condition1_concatenated.depthVert")

#input2 <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/subset/Results/Data/Intermediate-files/tiRNA.condition2_concatenated.depthVert")

df1 <- split( input1 , f = input1$V1 )  # Split dataframe based on column 1 elements
df2 <- split( input2 , f = input2$V1 )  # Split dataframe based on column 1 elements

results.df <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("feature",
                                                                 "Kolmogorov-smirnov.5prime.pvalue (cond1 vs cond2)",
                                                                 "Kolmogorov-smirnov.3prime.pvalue (cond1 vs cond2)",
                                                                 "-Log10.of.5prime.pvalue",
                                                                 "-Log10.of.3prime.pvalue",
                                                                 "Difference.of.-Log10")) # Initialise empty dataframe with column headers
count <- 1
for(subset1 in df1) {
  subset2 <- df2[[count]]  # Get the corresponding subset from the second file
  count <- count + 1
  feature <- as.character(subset1[1,1])
  subset1.length <- nrow(subset1)
  half.length <- as.integer(subset1.length/2)
  #if (half.length < 150) {
  #  length.penalty <- half.length/150###
  #} else {
  #  length.penalty <- 1
  #}
  ### Condition1
  cond1.fiveprime <- subset1[1:half.length,]
  cond1.five.distribution <- cond1.fiveprime$V3
  cond1.threeprime <- subset1[(half.length+1):subset1.length,]
  cond1.three.distribution <- cond1.threeprime$V3
  ### Condition2
  cond2.fiveprime <- subset2[1:half.length,]
  cond2.five.distribution <- cond2.fiveprime$V3
  cond2.threeprime <- subset1[(half.length+1):subset1.length,]
  cond2.three.distribution <- cond2.threeprime$V3
  ### KS tests to compare 5' and 3' in each condition. Get -Log10 of p-values
  ks.fiveprime <- ks.test(cond1.five.distribution, cond2.five.distribution)
  ks.fiveprime.pvalue <- ks.fiveprime$p.value
  ks.fiveprime.pvalue.negLog10 <- -log10(ks.fiveprime.pvalue) 
  ks.threeprime <- ks.test(cond1.three.distribution, cond2.three.distribution)
  ks.threeprime.pvalue <- ks.threeprime$p.value
  ks.threeprime.pvalue.negLog10 <- -log10(ks.threeprime.pvalue)
  if (ks.fiveprime.pvalue > 0.0) {  
    # Do nothing
  } else {
    ks.fiveprime.pvalue.negLog10 <- 20
  }
  if (ks.threeprime.pvalue > 0.0) {
    # Do nothing
  } else {
    ks.threeprime.pvalue.negLog10 <- 20
  }
  ### Compare difference
  maxval <- max(c(ks.fiveprime.pvalue.negLog10, ks.threeprime.pvalue.negLog10))
  minval <- min(c(ks.fiveprime.pvalue.negLog10, ks.threeprime.pvalue.negLog10))
  diff <- maxval - minval
  ### If the difference between the two is equivalent to a KS p-value of 0.001 (i.e. -log10 of 0.001 is 3), write to file
  #if ((diff) > 3) {
  results.df[nrow(results.df) + 1,] = list(feature, 
                                             ks.fiveprime.pvalue, 
                                             ks.threeprime.pvalue, 
                                             ks.fiveprime.pvalue.negLog10, 
                                             ks.threeprime.pvalue.negLog10,
                                             diff)
  #}
  

}

results.df <- results.df[order(-results.df$`Difference.of.-Log10`),]
newdata <- results.df[ which(results.df$`Difference.of.-Log10` > 3), ] # Equivalent to p-value of 0.001. More stringent than 0.05

write.table(results.df, 
            file = paste0(args[3], ".all-features.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
write.table(newdata, 
            file = paste0(args[3], ".potentially-cleaved-features.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

pdf.width <- nrow(newdata)*0.2 + 3
pdf(file = paste0(args[2], ".potentially-cleaved-features.pdf"), width = pdf.width, height = 5)
ggplot(data = newdata, mapping = aes(feature, newdata$`Difference.of.-Log10`, color=newdata$`Difference.of.-Log10`)) +
  geom_point() +
  ggtitle("Feature cleavage likelihood (-Log10 cutoff = 3)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_gradient(low="blue", high="red")
dev.off()


