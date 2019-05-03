#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script removes introns from genomecov file
###
###-------------------------------------------------------------------------------

### Use command line arguments
args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provided. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
}

### Read in files and get file basename
input1 <- read.table(args[1])
input2 <- read.table(args[2])

features <- input1$V1   # Group features by name
featuresUnion <- union(features, features)

for(feature in featuresUnion) {   # Loop over all tRNAs
  subset1 <- input1[grep(feature, input1$V1),] # Get the tRNA coordinates and count data for each tRNA
  if(any(input2==feature) == TRUE) {  # If the current tRNA is in input 2 fie, it has an intron and this must be removed
    myRow <- input2[grep(feature, input2$V1),]  # Get the row that contains this intron-containing tRNA information
    trna <- as.character(myRow$V1)  ##
    start <- as.numeric(myRow$V2)   ## Get the intron information
    stop <- as.numeric(myRow$V3)    ##
    newdf <- subset1[-c(start:stop),]   # Remove the rows that correspond to the intron
    newdf$V2 <- 1:nrow(newdf)    #Replace column 2 with range 1 to dataframe end
    write.table(newdf, 
                file = args[3], 
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE,
                sep = "\t",
                append = TRUE)
  } else {
    write.table(subset1, 
                file = args[3], 
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE,
                sep = "\t",
                append = TRUE)
    }
  }
