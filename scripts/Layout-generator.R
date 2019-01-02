#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  message("No arguments provided. Defaulting to condition 1 vs condition 2 assuming three replicates each.")
} else if (length(args)>2) {
  stop("Error: Too many command line arguments. Quitting.")
} else if (length(args)==1) { 
  stop("Error: Require 2 command line arguments e.g. 'Layout-generator.R /path/to/FASTQ-files/ /output/directory/' ")
} else if (length(args)==2) { 
  myPath <- args[1]
  outDir <- args[2]
  if (dir.exists(myPath)) {
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
    setwd(outDir)
    write.table(df2, sep = ",", quote = FALSE, file = "predicted_exp_layout.csv", row.names = FALSE)
    message("\nClustering resulted in this division (see column 3):")
    rownames(df2) <- NULL
    print(df2)
    message(paste0("\nIf these groupings are incorrect, please rearrange the 'predicted_exp_layout.csv' file\nso that replicates are grouped together.\n\nWriting output file to ", outDir))
  } else {
    stop("Error: Command line argument 1 is not a directory path or this directory does not exist.")
  }
}