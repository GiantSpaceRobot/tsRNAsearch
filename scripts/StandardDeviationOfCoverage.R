#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script calculates standard deviation of percentage of read coverage difference
### 
###-------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
} 

#library(data.table)

#if (length(args)==4) {
#  GTF <- read.table(args[4], sep = "\t")
#GTF <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/DBs/hg19-snomiRNA.gtf", sep = "\t")
#} 

### https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#### Input file base names
name1 <- "condition1"
name2 <- "condition2"
#name1 <- basename(args[1])
#name2 <- basename(args[2])

#### Input file
input <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/TOPHAT-FullTest/Results/Data/stddev.txt")
input <- read.table(args[1])

plot(input$V4)

features <- input$V1   # Group features by name
featuresUnion <- union(features, features) # Get unique set of features

#pdf(args[2])
#pdf("/home/paul/Documents/Pipelines/tirna-pipeline/TOPHAT-FullTest/Results/Plots/CoverageVariability.pdf")

my.list = list("Feature", "Mean of mean1-mean2", "Std.dev of mean1-mean2", "Coef. of var of mean1-mean2", "Skewness of mean1-mean2", "Fragment retention higher in:")
write.table(my.list, 
            file = "/home/paul/Documents/Pipelines/tirna-pipeline/TOPHAT-FullTest/Results/Data/Information-2.txt",
            #file = args[2],
            quote = FALSE, 
            append = TRUE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

for(feature in featuresUnion) {
  #feature <- "AsnATT"
  #subset <- input[grep("AsnATT", input$V1),]
  subset <- input[grep(feature, input$V1),]
  subset$V5 <- (subset$V2-subset$V3)
  mean.of.difference <- mean(subset$V5)
  #geometric.mean <- gm_mean(subset$V5)
  geometric.mean <- exp(mean(log(subset$V5)))
  geometric.sd <- exp(sd(log(subset$V5)))
  geometric.coef.of.var <- geometric.sd**(1/geometric.mean)   # Geometric coef.of.var = geo.sd^(1/geo.mean)
  #stddev.of.percent <- sd(subset$V4)
  stddev.of.difference <- sd(subset$V5)
  coef.of.variation <- (stddev.of.difference/mean.of.difference * 100)
  median.of.difference <- median(subset$V5)
  skewness <- 3*(median.of.difference/mean.of.difference)/stddev.of.difference
  if (mean.of.difference > 0) {
    condition <- 1
  } else {
    condition <- 2
  }
  my.list <- list(feature, mean.of.difference, stddev.of.difference, coef.of.variation, skewness, paste0("Condition ", condition))
  write.table(my.list, 
              file = "/home/paul/Documents/Pipelines/tirna-pipeline/TOPHAT-FullTest/Results/Data/Information-2.txt",
              #file = args[2],
              append = TRUE, 
              quote = FALSE, 
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  #plot(subset$V4, 
  #     type = "l",
  #     col = "blue", 
  #     lwd = 2, 
  #     main = feature,
  #     xlab = "Nucleotide position",
  #     ylab = "Coverage variability (%)")
}
#dev.off()


