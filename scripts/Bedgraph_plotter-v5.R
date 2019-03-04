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

if (length(args)==4) {
  GTF <- read.table(args[4], sep = "\t")
  #GTF <- read.table("/home/paul/Documents/Pipelines/tsRNAsearch/DBs/hg19-snomiRNA_cdhit.gtf", sep = "\t")
} 

### Change this if you want only to plot features with a higher mean coverage
### (e.g. mean coverage of 100 reads is 'mean.cutoff <- 100')
mean.cutoff <- as.integer(args[3])  
#mean.cutoff <- 20

### Input file
fileData <- read.table(args[1], 
                       col.names = c("chrom", "coordinate", "dataValue"))
#fileData <- read.table("/home/paul/Documents/Pipelines/tsRNAsearch/test/tRNA-alignment/accepted_hits_sorted.depth", col.names = c("chrom", "coordinate", "dataValue"))
#feature <- "ENSG00000207730"
#feature <- "ValCAC"

#features <- fileData$chrom   # Group features by name
#featuresUnion <- union(features, features) # Get unique set of features

df1 <- split( fileData , f = fileData$chrom )  # Split dataframe based on column 1 elements


pdf(args[2])
#pdf("/home/paul/Documents/Pipelines/tirna-pipeline/output.pdf3")
for(subset in df1) {
  feature <- as.character(subset[1,1])
  #subset.NoFlanks <- tail(subset, -10)
  #subset.NoFlanks <- head(subset.NoFlanks, -10)
  subset.NoFlanks <- subset
  subset.mean <- mean(subset.NoFlanks$dataValue)
  if (length(args)==4) {
    featureRows <- GTF[grep(feature, GTF$V9),]
    featureRows <- featureRows[1,]
    geneName <- as.character(sub(".*gene_name *(.*?) *; gene_source.*", "\\1", featureRows$V9))
    feature <- geneName
  } 
  if(subset.mean>mean.cutoff){
    plot(subset.NoFlanks$dataValue, 
       type = "l",
       col = "blue", 
       lwd = 2, 
       main = feature,
       xlab = "Nucleotide position",
       ylab = "Coverage depth (Reads Per Million)")
  }
}
dev.off()


