#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Require 2 command line arguments e.g. 'Rscript Layout-generator.R /path/to/FASTQ-files/ Layout.csv' ")
} else if (length(args)>2) {
  stop("Require 2 command line arguments e.g. 'Rscript Layout-generator.R /path/to/FASTQ-files/ Layout.csv' ")
} else if (length(args)==1) { 
  stop("Require 2 command line arguments e.g. 'Rscript Layout-generator.R /path/to/FASTQ-files/ Layout.csv' ")
} else if (length(args)==2) { 
  myPath <- args[1]
  #outDir <- args[2]
  ### Make sure path name starts and ends in slash (need full path name)
  #if (startsWith(args[2], "/") == TRUE) {
  #} else {
  #stop("Output directory path must be full path and must begin with '/'")
  #}

  if (dir.exists(myPath)) {
    outDir <- getwd()
    setwd(myPath)
    message(paste0("\nCarrying out hierarchical clustering of filenames to classify files into groups."))
    ### A weird idea of clustering the filenames based on string similarity 
    ### (so I don't need to ask the user which samples belong to the same group)
    file.names <- dir(pattern=".")
    d <- adist(file.names)
    rownames(d) <- file.names
    hc <- hclust(as.dist(d))
    df <- data.frame(file.names,cutree(hc,k=2))
    df2 <- df
    colnames(df2) <- NULL
    df2[,2] <- paste0("Condition",df2[,2],sep="")
    #setwd(outDir)
    #write.table(df2, sep = ",", quote = FALSE, file = "predicted_exp_layout.csv", row.names = FALSE)
    write.table(df2, sep = ",", quote = FALSE, file = paste0(outDir, "/", args[2], sep=""), row.names = FALSE)
    message("\nClustering resulted in this division (see column 3):")
    rownames(df2) <- NULL
    print(df2)
    message(sprintf("\nIf these groupings are incorrect, please rearrange the '%s' file\nso that replicates are grouped together.\n", args[2]))
  } else {
    stop("Command line argument 1 is not a directory path or this directory does not exist.")
  }
}
