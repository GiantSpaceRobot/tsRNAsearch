#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script calculates mean and stdev of read coverage files (.genomecov)
### 
###-------------------------------------------------------------------------------

### Use command line arguments
args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provided. If not, return an error.
if (length(args)==0) {
  stop("Not enough command line arguments provided. Input file and output file names required.")
}

### Read in files and get file basename
input1 <- read.table(args[1])
### Get every third column to calculate mean and stdev (this loop can deal with variable numbers of replicates)
input1$mean <- rowMeans(subset(input1, select = c(rep(FALSE,2), TRUE)), na.rm = TRUE)
input1$stdev <- apply(input1[c(rep(FALSE,2), TRUE)], 1, sd) # Calculate stdev using every third column
### Write file
write.table(x = input1,
            sep = "\t",
            file = args[2], 
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE)