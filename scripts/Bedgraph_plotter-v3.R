#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script makes coverage plots
### 
###-------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
}

### Change this if you want only to plot features with a higher mean coverage
### (e.g. mean coverage of 100 reads is 'mean.cutoff <- 100')
mean.cutoff <- 0   

### Input file
fileData <- read.table(args[1], 
                       col.names = c("chrom", "coordinate", "dataValue"))

features <- fileData$chrom   # Group features by name
featuresUnion <- union(features, features) # Get unique set of features

pdf(args[2])
for(feature in featuresUnion) {
  subset <- fileData[grep(feature, fileData$chrom),]
  subset.NoFlanks <- tail(subset, -10)
  subset.NoFlanks <- head(subset.NoFlanks, -10)
  subset.mean <- mean(subset.NoFlanks$dataValue)
  if(subset.mean>mean.cutoff){
    plot(subset.NoFlanks$dataValue, 
       type = "l",
       col = "blue", 
       lwd = 2, 
       main = feature,
       xlab = "Nucleotide position",
       ylab = "Coverage depth")
  }
}
dev.off()


