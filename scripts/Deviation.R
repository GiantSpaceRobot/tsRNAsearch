#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script calculates standard deviation of read coverage across features and compares conditions.
### 
###-------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
} 

library(data.table)

#if (length(args)==4) {
#  GTF <- read.table(args[4], sep = "\t")
#GTF <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/DBs/hg19-snomiRNA.gtf", sep = "\t")
#} 

#### Input file base names
name1 <- "condition1"
name2 <- "condition2"
name1 <- basename(args[1])
name2 <- basename(args[2])

#### Input file
input1 <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/TOPHAT-FullTest/Results/Data/condition1_concatenated.tiRNA.genomecov")
input2 <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/TOPHAT-FullTest/Results/Data/condition2_concatenated.tiRNA.genomecov")
input1 <- read.table(args[1])
input2 <- read.table(args[2])

### Create header for file input1
col.num1 <- as.integer(length(input1))
header1 <- vector("list", col.num1)
for(i in seq(1, col.num1)) {
  header1[[i]] <- paste0(name1, ".V", i)
}

### Create header for file input2
col.num2 <- as.integer(length(input2))
header2 <- vector("list", col.num2)
for(i in seq(1, col.num2)) {
  header2[[i]] <- paste0(name2, ".V", i)
}

### Add headers to dataframes
colnames(input1) <- header1
colnames(input2) <- header2

### Get sum, mean and std of data
sum1.colname <- paste0(name1, ".sum")
mean1.colname <- paste0(name1, ".mean")
stddev1.colname <- paste0(name1, ".stddev")


for (row in 1:nrow(input1)) {
  for (columns in 1:ncol(row)) {
    if(columns%%3) {   # If column is 3rd, 6th, 9th etc. (divisible by 3)
      print(paste(row, " Column is: ", columns))
    }
  }
}

#input1$sum1.colname <- colSums()


### Merge dataframes side-by-side
new.df <- cbind(input1, input2)




#                       col.names = c("chrom", "coordinate", "dataValue"))
##fileData <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/SimReads_HISAT/snomiRNA-alignment/accepted_hits.genomecov", col.names = c("chrom", "coordinate", "dataValue"))
##feature <- "ENSG00000207730"
#
#features <- fileData$chrom   # Group features by name
#featuresUnion <- union(features, features) # Get unique set of features
#
#pdf(args[2])
##pdf("/home/paul/Documents/Pipelines/tirna-pipeline/output.pdf3")
#for(feature in featuresUnion) {
#  subset <- fileData[grep(feature, fileData$chrom),]
#  subset.NoFlanks <- tail(subset, -10)
#  subset.NoFlanks <- head(subset.NoFlanks, -10)
#  subset.mean <- mean(subset.NoFlanks$dataValue)
#  if (length(args)==4) {
#    featureRows <- GTF[grep(feature, GTF$V9),]
#    featureRows <- featureRows[1,]
#    geneName <- as.character(sub(".*gene_name *(.*?) *; gene_source.*", "\\1", featureRows$V9))
#    feature <- geneName
#  } 
#  if(subset.mean>mean.cutoff){
#    plot(subset.NoFlanks$dataValue, 
#         type = "l",
#         col = "blue", 
#         lwd = 2, 
#         main = feature,
#         xlab = "Nucleotide position",
#         ylab = "Coverage depth")
#  }
#}
#dev.off()


