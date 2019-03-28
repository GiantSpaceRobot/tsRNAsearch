#!/usr/bin/env Rscript

## Script name: DESeq2_tiRNA-pipeline.R
##
## Purpose of script: Carry out differential gene expression analysis of sncRNAs/tiRNAs from small RNA-seq data
##
## Author: Dr. Paul Donovan
##
## Date Created: 2-1-2018 (2nd of Jan 2018)
##
## Version: 10
##
## Copyright: MIT license
##
## Email: pauldonovan@rcsi.com

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  message("No arguments provided. Defaulting to condition 1 vs condition 2 assuming three replicates each.")
} else if (length(args)>2) {
  stop("Error: Too many command line arguments. Quitting.")
} else if (length(args)==1) {             # If only a filepath was provided:
  myPath <- args[1]
  myPath <- "/home/paul/Documents/Pipelines/tirna-pipeline/Full-test/Results/Data/Intermediate-files/Counts/"
  if (dir.exists(myPath)) {
    setwd(myPath)
    message("Carrying out hierarchical clustering of filenames to classify files into groups.")
    ### A weird idea of clustering the filenames based on string similarity 
    ### (so I don't need to ask the user which samples belong to the same group)
    file.names.long <- dir(pattern=".all-features.count")
    file.names <- gsub(".all-features.count", "", file.names.long)
    d <- adist(file.names)
    rownames(d) <- file.names
    hc <- hclust(as.dist(d))
    ### Create directory for the results files
    if (file.exists("DE_Results")){
    } else {
      dir.create("DE_Results")
    }
    pdf("DE_Results/Hierarchical-Clustering-of-Filenames.pdf",
        width=12,height=12)
    plot(hc, cex=.4)
    rect.hclust(hc,k=2)
    text(2,10, "If these groupings are incorrect,\nplease provide a CSV file indicating\nthe correct groups (see README)\nand restart the analysis")
    dev.off()
    df <- data.frame(file.names,cutree(hc,k=2))
    df2 <- df
    colnames(df2) <- NULL
    df2[,2] <- paste0("Condition",df2[,2],sep="")
    write.table(df2, sep = ",", quote = FALSE, file = "predicted_exp_layout.csv", row.names = FALSE)
    lvls.df <- as.data.frame(table(df$cutree.hc..k...2.))
    Condition1 <- "Condition1"
    ReplicateNumber1 <- lvls.df[1,2]
    Condition2 <- "Condition2"
    ReplicateNumber2 <- lvls.df[2,2]
  } else {
    stop("Error: Command line argument 1 is not a directory path.")
  }
} else if (length(args)==2) {     # If a .csv and file path were provided:
  experiment <- args[1]
  exp.file <- read.csv(experiment, header=FALSE)
  setwd(args[2])
  
  ### Count the number of factor levels in the provided .csv file (this should be 2)
  lvls <- levels(exp.file$V2)
  if (length(lvls)==2) {
    #Do nothing
  } else {
    lvls.pasted <- paste(unlist(lvls), collapse=',')
    quit.message <- sprintf("Wrong number of factor levels: %s. There should be two factor levels", lvls.pasted)
    stop(quit.message)
  }
  ### Create directory for the results files
  if (file.exists("DE_Results")){
  } else {
    dir.create("DE_Results")
  }
  file.CSV <- read.csv(args[1], header=FALSE)
  lvls.df <- as.data.frame(table(file.CSV$V2))
  ReplicateNumber1 <- lvls.df[1,2]
  Condition1 <- toString(lvls.df[1,1])
  ReplicateNumber2 <- lvls.df[2,2]
  Condition2 <- toString(lvls.df[2,1])
  myPath <- args[2]
  message("A .csv file was provided, and the directory provided exists.")
}

### Make sure path name starts and ends in slash (need full path name)
if (startsWith(myPath, "/") == TRUE) {
} else {
  stop("Directory path must be full path and must begin with '/'")
}
### Must end in slash
if (substring(myPath, nchar(myPath)) == "/") {
} else {
  myPath <- paste0(myPath, "/")
}

### A variable that will act as a prefix to all output files
ResultsFile <- paste0(Condition1,"-vs-",Condition2)


#-----------------#
#    Libraries    #
#-----------------#

# If libraries not installed, install them
#source("http://bioconductor.org/biocLite.R")
#biocLite(DESeq2)

library(DESeq2)
library(gplots)
library(ggplot2)


#------------------------------#
#    Define DESeq2 function    #
#------------------------------#

### Function to carry out DESeq2 analysis on provided files
DESeq2.function <- function(path.to.files){
  
  #path.to.files <- myPath
  ### Create checkpoints (for simple degbugging) 
  count <- 0
  count <- count + 1
  print(paste0("Checkpoint ", count))
  
  ### Read files
  file.names <- dir(path.to.files, pattern =".count")
  cDataAll <- NULL
  for (i in 1:length(file.names)){
    full.path <- paste0(path.to.files, file.names[i])
    file <- read.table(full.path, header = TRUE)
    cDataAll <- cbind(cDataAll, file[,2])
  }
  colnames(cDataAll) <- (file.names)
  rownames(cDataAll) <- file[,1]
  
  ### checkpoint
  count <- count + 1
  print(paste0("Checkpoint ", count))
  
  ###
  groups <- factor(x=c(rep(Condition1, ReplicateNumber1), rep(Condition2, ReplicateNumber2)), levels=c(Condition1, Condition2))
  tpm <- t(t(cDataAll)/colSums(cDataAll))*1e6
  inlog <- log(tpm)
  colLabel <- c(rep("#E41A1C", ReplicateNumber1), rep("#377EB8", ReplicateNumber2))
  colTy <- c(rep(1:ReplicateNumber1, ReplicateNumber1), rep(1:ReplicateNumber2, ReplicateNumber2))
  
  ### checkpoint
  count <- count + 1
  print(paste0("Checkpoint ", count))

  ### Create a density plot
  #####pdf(paste0("DE_Results/", ResultsFile, "_Density-Plot.pdf"),
  #####    width=12,height=12)
  #####plot(density(inlog[,1]), 
  #####     ylim=c(0,0.4), 
  #####     main="Density plot of counts per gene", 
  #####     lty=colTy[1], 
  #####     xlab="Log of TPM per gene", 
  #####     ylab="Density", 
  #####     col=colLabel[1])
  #####for(i in 2:ncol(tpm)){
  ##### lines(density(inlog[,i]), lty=colTy[i], col=colLabel[i])
  #####}
  #####legend("topright", legend=colnames(tpm), lty=colTy, col=colLabel)
  #####dev.off()
  
  ### checkpoint
  count <- count + 1
  print(paste0("Checkpoint ", count))

  ### Differential Expression Analysis
  colData <- DataFrame(condition=groups) 
  dds <- DESeqDataSetFromMatrix(cDataAll, colData, formula(~condition))
  dds <- DESeq(dds)
  res <- results(dds, cooksCutoff=FALSE)
  res <- res[order(res$log2FoldChange, decreasing=TRUE),]
  
  ### checkpoint
  count <- count + 1
  print(paste0("Checkpoint ", count))

  ### Carry out a log transformation of the DESeq2 data set
  rld <- rlogTransformation(dds)

  ### Create a histogram
  #####pdf(paste0("DE_Results/", ResultsFile, "_Histogram.pdf"),
  #####    width=12,height=12)
  #####hist(assay(rld))
  #####dev.off()
  
  ### checkpoint
  count <- count + 1
  print(paste0("Checkpoint ", count))

  ### Convert the rld transformation into a matrix
  rld.matrx <- assay(rld)   
  rld.df <- data.frame(rld.matrx)
  
  ### Sample distance heatmap
  sampleDists <- as.matrix(dist(t(assay(rld))))
  #(mycols <- brewer.pal(length(file.names), "Paired")[1:length(unique(file.names))])
  #####pdf(paste0("DE_Results/", ResultsFile, "_Distance-Matrix.pdf"),
  #####    width=12,height=12)
  #####heatmap.2(as.matrix(sampleDists), key=F, trace="none",
  #####          col=colorpanel(100, "black", "white"),
  #####          #ColSideColors=mycols[file.names], 
  #####          #RowSideColors=mycols[file.names],
  #####          cexRow = 0.8,
  #####          cexCol = 0.8,
  #####          margin=c(10, 13), 
  #####          srtCol=45,
  #####          main="Sample Distance Matrix")
  #####dev.off()
  
  ### Feature heatmap
  topVarGenes <- head(order(rowVars(assay(rld)), decreasing=T),50)
  matrix.heatmap <- assay(rld)[ topVarGenes, ]
  matrix.heatmap <- matrix.heatmap - rowMeans(matrix.heatmap)
  heatmap.2(matrix.heatmap)
  # select the 'contrast' you want
  annotation_data <- as.data.frame(colData(rld)[c(Condition1,Condition2)])
  pheatmap(matrix.heatmap, annotation_col=annotation_data)
  
  ### checkpoint
  count <- count + 1
  print(paste0("Checkpoint ", count))
  
  ### Principal component analysis
  #Tpm PCA plot (not DESeq2)
  d <- dist(t(tpm))
  fit=cmdscale(d, eig=TRUE, k=2)
  x=fit$points[,1]
  y=fit$points[,2]
  #####pdf(paste0("DE_Results/", ResultsFile, "_tpm-PCA.pdf"),
  #####    width=8,height=8)
  #####par(xpd = T, mar = par()$mar + c(5,4,4,8))
  #####plot(x, y, 
  #####     type="p", 
  #####     pch=20, 
  #####     col=colLabel,
  #####     bty="L")
  #####box()
  #####legend("right", 
  #####       inset=c(-0.7,0), 
  #####       pch=20,
  #####       col=colLabel,
  #####       cex = 0.7, 
  #####       legend=colnames(tpm),
  #####       xpd = TRUE,
  #####       bty = "n",
  #####       lty=NULL)
  #####dev.off()
  #par(mar=c(5, 4, 4, 2) + 0.1)
  
  ### checkpoint
  count <- count + 1
  print(paste0("Checkpoint ", count))

  # PCA plot that won't print to PDF
  #pdf(paste0("DE_Results/", ResultsFile, "_PCA.pdf"),
  #    width=12,height=12)
  #pcaData <- plotPCA(rld, returnData=TRUE)
  #percentVar <- round(100 * attr(pcaData, "percentVar"))
  #ggplot(pcaData, aes(PC1, PC2, colour=name)) +
  #  geom_point(size=3) +
  #  geom_text(aes(label=group),hjust=-0.25, vjust=0.2) +
  #  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #  coord_fixed()
  #dev.off()
  
  ### Write DESeq2 results to file
  #####if (file.exists("DE_Results/DESeq2")){
  #####} else {
  #####  dir.create("DE_Results/DESeq2")
  #####}
  #####
  ######## checkpoint
  #####count <- count + 1
  #####print(paste0("Checkpoint ", count))
#####
  #####write.csv(as.data.frame(res), file=paste0("DE_Results/DESeq2/", ResultsFile, "_DESeq2-output.csv"))
  #####
  #####up <- (res[!is.na(res$padj) & res$padj <= 0.1 &
  #####                    res$log2FoldChange >= 0.5, ])
  #####
  #####down <- (res[!is.na(res$padj) & res$padj <= 0.1 &    
  #####                      res$log2FoldChange <= -0.5, ])    
  ######print(sprintf("%s genes up-regulated, %s genes down-regulated", length(up), length(down)))
  #####write.table(x = up,
  #####            sep = ",",
  #####            file=paste0("DE_Results/DESeq2/", ResultsFile, "_DESeq2-output-upregulated.csv"), 
  #####            col.names=NA,
  #####            row.names=TRUE, 
  #####            quote = FALSE)
  #####write.table(x = down,
  #####            sep = ",",
  #####            file=paste0("DE_Results/DESeq2/", ResultsFile, "_DESeq2-output-downregulated.csv"), 
  #####            col.names=NA,
  #####            row.names=TRUE, 
  #####            quote = FALSE)
  #####
  ### checkpoint
  count <- count + 1
  print(paste0("Checkpoint ", count))

  }


#--------------------#
#    Run analysis    #
#--------------------#

DESeq2.function(myPath)

