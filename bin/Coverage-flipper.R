#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script flips the read counts of features in the minus orientation (genomecov file + gtf for orientation info)
### i.e. if a feature is on the minus strand of DNA, make it appear to be on the plus strand for collapsing into one tRNA
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

for(feature in featuresUnion) {
  i <- grepl(paste0("^", feature), input2$V1)   # Get row from input2 containing the tRNA
  newdf <- input2[i,]
  if(newdf$V7 == "-") {   # If this tRNA is on the minus strand, flip read counts
    minus.tRNAs <- grepl(paste0("^", feature), input1$V1)   # Get row from input2 containing the tRNA
    minus.tRNAs.df <- input1[minus.tRNAs,]
    minus.tRNAs.df$V3 <- rev(minus.tRNAs.df$V3)  # Reverse column 3 so the data reflects plus strand data
    write.table(minus.tRNAs.df, 
                file = args[3], 
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE,
                sep = "\t",
                append = TRUE)
  } else {    # This tRNA is on the plus strand, write to file as is
    plus.tRNAs <- grepl(paste0("^", feature), input1$V1)   # Get row from input2 containing the tRNA
    plus.tRNAs.df <- input1[plus.tRNAs,]
    write.table(plus.tRNAs.df, 
                file = args[3], 
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE,
                sep = "\t",
                append = TRUE)
  }
}
