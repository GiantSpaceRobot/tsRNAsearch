#!/usr/bin/env Rscript

## Script name: DESeq2_tsRNAsearch.R
##
## Purpose of script: Carry out differential gene expression analysis of ncRNAs/genes from small RNA-seq data
##
## Author: Dr Paul Donovan
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
} else if (length(args)>3) {
  stop("Error: Too many command line arguments. Quitting.")
} else if (length(args)==1) {
  stop("Error: Require two command line arguments at least: path-to-files, output-prefix, [layout-file.csv]. Quitting.")
} else if (length(args)==2) {             # If only a filepath was provided:
  myPath <- args[1]
  if (dir.exists(myPath)) {
    message("Carrying out hierarchical clustering of filenames to classify files into groups.")
    ### A weird idea of clustering the filenames based on string similarity 
    ### (so I don't need to ask the user which samples belong to the same group)
    file.names.long <- dir(path = myPath, pattern=".all-features.count")
    file.names <- gsub(".all-features.count", "", file.names.long)
    d <- adist(file.names)
    rownames(d) <- file.names
    hc <- hclust(as.dist(d))
    ### Create directory for the results files
    if (file.exists("DE_Results")){
    } else {
      dir.create("DE_Results")
    }
    pdf(paste0(myPath, "DE_Results/Hierarchical-Clustering-of-Filenames.pdf"),
        width=12,height=12)
    plot(hc, cex=.4)
    rect.hclust(hc,k=2)
    text(2,10, "If these groupings are incorrect,\nplease provide a CSV file indicating\nthe correct groups (see README)\nand restart the analysis")
    garbage <- dev.off()
    df <- data.frame(file.names,cutree(hc,k=2))
    df2 <- df
    colnames(df2) <- NULL
    df2[,2] <- paste0("Condition",df2[,2],sep="")
    write.table(df2, sep = ",", quote = FALSE, file = paste0(myPath, "predicted_exp_layout.csv"), row.names = FALSE)
    lvls.df <- as.data.frame(table(df$cutree.hc..k...2.))
    Condition1 <- "Condition1"
    ReplicateNumber1 <- lvls.df[1,2]
    Condition2 <- "Condition2"
    ReplicateNumber2 <- lvls.df[2,2]
  } else {
    stop("Error: Command line argument 1 is not a directory path.")
  }
} else if (length(args)==3) {     # If a .csv and file path were provided:
  experiment <- args[3]
  exp.file <- read.csv(experiment, header=FALSE)
  myPath <- args[1]

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
  if (file.exists(paste0(myPath, "DE_Results"))){
  } else {
    dir.create(paste0(myPath, "DE_Results"))
  }
  file.CSV <- read.csv(args[3], header=FALSE)
  lvls.df <- as.data.frame(table(file.CSV$V2))
  ReplicateNumber1 <- lvls.df[1,2]
  Condition1 <- toString(lvls.df[1,1])
  ReplicateNumber2 <- lvls.df[2,2]
  Condition2 <- toString(lvls.df[2,1])
  myPath <- args[1]
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

### Directory for writing DESeq2 results to file
if (file.exists(paste0(myPath, "DE_Results/DESeq2"))){
} else {
  dir.create(paste0(myPath, "DE_Results/DESeq2"))
}

### If replicate number = 1, stop analysis
if(ReplicateNumber1==1) {
  file.create(paste0(myPath, "DE_Results/DESeq2/upregulated.csv"))
  file.create(paste0(myPath, "DE_Results/DESeq2/downregulated.csv"))
  quit.message <- "Replicate number is 1, cannot continue DESeq2 analysis"
  stop(quit.message)
}

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
  
  ### Create checkpoints (for simple degbugging) 
  print("Checkpoint 1")
  
  ### Read files
  path.to.files <- myPath
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
  print("Checkpoint 2")
  
  ###
  groups <- factor(x=c(rep(Condition1, ReplicateNumber1), rep(Condition2, ReplicateNumber2)), levels=c(Condition1, Condition2))
  tpm <- t(t(cDataAll)/colSums(cDataAll))*1e6
  inlog <- log(tpm)
  colLabel <- c(rep("#E41A1C", ReplicateNumber1), rep("#377EB8", ReplicateNumber2))
  colTy <- c(rep(1:ReplicateNumber1, ReplicateNumber1), rep(1:ReplicateNumber2, ReplicateNumber2))
  
  ### checkpoint
  print("Checkpoint 3")

  ### Create a density plot
  pdf(paste0(path.to.files, "DE_Results/", ResultsFile, "_Density-Plot.pdf"),
      width=12,height=12)
  plot(density(inlog[,1]), 
       ylim=c(0,0.4), 
       main="Density plot of counts per gene", 
       lty=colTy[1], 
       xlab="Log of TPM per gene", 
       ylab="Density", 
       col=colLabel[1])
  for(i in 2:ncol(tpm)){
   lines(density(inlog[,i]), lty=colTy[i], col=colLabel[i])
  }
  legend("topright", legend=colnames(tpm), lty=colTy, col=colLabel)
  garbage <- dev.off()
  
  ### checkpoint
  print("Checkpoint 4")

  ### Differential Expression Analysis
  colData <- DataFrame(condition=groups) 
  dds <- DESeqDataSetFromMatrix(cDataAll, colData, formula(~condition))
  dds <- DESeq(dds)
  res <- results(dds, cooksCutoff=FALSE)
  res <- res[order(res$log2FoldChange, decreasing=TRUE),]
  
  ### checkpoint
  print("Checkpoint 5")

  ### Carry out a log transformation of the DESeq2 data set
  tryCatch({
    print("Attempting variance stabilising transformation of dds")
    #I hit an error with simulated data and rlog transformation. Michael Love (author) recommended VST https://support.bioconductor.org/p/100927/
    rld <- varianceStabilizingTransformation(dds)
  }, warning = function(w) {
    # Do nothing
  }, error = function(e) {
    print("Failed variance stabilising transformation. Using rlog transformation instead")
    rld <- rlogTransformation(dds)
  }, finally = {
    # Do nothing
  })
  if(exists("rld")==TRUE){
    #Do nothing
  } else {
    print("
          
          WARNING: COULD NOT TRANSFORM COUNT DATA AS IT IS IRREGULAR. CONTINUING DESEQ2 ANALYSIS USING UNTRANSFORMED DATA.
          
          ")
    rld <- dds
    #stop("Error: Could not apply transformation to data. Aborted DESeq2 analysis.")
  }
  
  ### Create a histogram
  pdf(paste0(path.to.files, "DE_Results/", ResultsFile, "_Histogram.pdf"),
      width=12,height=12)
  hist(assay(rld))
  garbage <- dev.off()
  
  ### checkpoint
  print("Checkpoint 6")

  ### Convert the rld transformation into a matrix
  rld.matrx <- assay(rld)   
  rld.df <- data.frame(rld.matrx)
  
  ### Sample distance heatmap
  sampleDists <- as.matrix(dist(t(assay(rld))))
  #(mycols <- brewer.pal(length(file.names), "Paired")[1:length(unique(file.names))])
  pdf(paste0(path.to.files, "DE_Results/", ResultsFile, "_Distance-Matrix.pdf"),
      width=12,height=12)
  heatmap.2(as.matrix(sampleDists), key=F, trace="none",
            col=colorpanel(100, "black", "white"),
            #ColSideColors=mycols[file.names], 
            #RowSideColors=mycols[file.names],
            cexRow = 0.8,
            cexCol = 0.8,
            margin=c(10, 13), 
            srtCol=45,
            main="Sample Distance Matrix")
  garbage <- dev.off()
  
  ### checkpoint
  print("Checkpoint 7")
  
  ### Principal component analysis
  #Tpm PCA plot (not DESeq2)
  d <- dist(t(tpm))
  fit=cmdscale(d, eig=TRUE, k=2)
  x=fit$points[,1]
  y=fit$points[,2]
  pdf(paste0(path.to.files, "DE_Results/", ResultsFile, "_tpm-PCA.pdf"),
      width=8,height=8)
  par(xpd = T, mar = par()$mar + c(5,4,4,8))
  plot(x, y, 
       type="p", 
       pch=20, 
       col=colLabel,
       bty="L")
  box()
  legend("right", 
         inset=c(-0.7,0), 
         pch=20,
         col=colLabel,
         cex = 0.7, 
         legend=colnames(tpm),
         xpd = TRUE,
         bty = "n",
         lty=NULL)
  garbage <- dev.off()
  #par(mar=c(5, 4, 4, 2) + 0.1)
  
  ### checkpoint
  print("Checkpoint 8")

  write.csv(as.data.frame(res), file=paste0(path.to.files, "DE_Results/DESeq2/", ResultsFile, "_DESeq2-output.csv"))
  
  up <- (res[!is.na(res$padj) & res$padj <= 0.1 &
                      res$log2FoldChange >= 0.5, ])
  
  down <- (res[!is.na(res$padj) & res$padj <= 0.1 &    
                        res$log2FoldChange <= -0.5, ])    
  write.table(x = up,
              sep = ",",
              file=paste0(path.to.files, "DE_Results/DESeq2/", ResultsFile, "_DESeq2-output-upregulated.csv"), 
              col.names=NA,
              row.names=TRUE, 
              quote = FALSE)
  write.table(x = down,
              sep = ",",
              file=paste0(path.to.files, "DE_Results/DESeq2/", ResultsFile, "_DESeq2-output-downregulated.csv"), 
              col.names=NA,
              row.names=TRUE, 
              quote = FALSE)
  
  ### checkpoint
  print("Checkpoint 9")
  
  up.df <- as.data.frame(up)
  down.df <- as.data.frame(down)
  newdata.subset <- rbind(up.df, down.df)
  newdata.subset$features <- rownames(newdata.subset)
  newdata.subset$negLog10 <- -log10(newdata.subset$padj)
  tsRNAs.df <- subset(newdata.subset, startsWith(as.character(features), "chr"))
  # If there are more than 20 features, show top 20
  if(nrow(tsRNAs.df) > 20){
    tsRNAs.df.subset <- head(tsRNAs.df, n = 20)
  } else {
    tsRNAs.df.subset <- tsRNAs.df
  }
  genes.df <- subset(newdata.subset, startsWith(as.character(features), "ENS"))
  # If there are more than 20 features, show top 20
  if(nrow(genes.df) > 20){
    genes.df.subset <- head(genes.df, n = 20)
  } else {
    genes.df.subset <- genes.df
  }
  ### replace Inf values (extremely low adj p-value) with -Log10 of 300
  tsRNAs.df.subset <- data.frame(lapply(tsRNAs.df.subset, function(x) {gsub(Inf, "300", x)}))
  tsRNAs.df.subset$negLog10 <- as.numeric(levels(tsRNAs.df.subset$negLog10))[tsRNAs.df.subset$negLog10] #Convert factor type to numeric
  genes.df.subset <- data.frame(lapply(genes.df.subset, function(x) {gsub(Inf, "300", x)}))
  genes.df.subset$negLog10 <- as.numeric(levels(genes.df.subset$negLog10))[genes.df.subset$negLog10] #Convert factor type to numeric
  
  ### checkpoint
  print("Checkpoint 10")
  
  pdf.width <- nrow(tsRNAs.df.subset)*0.2 + 3
  pdf(file = paste0(path.to.files, args[2], "_tsRNAs.high-DE-negLog10_padj.pdf"), width = pdf.width, height = 5)
  print(ggplot(data = tsRNAs.df.subset, mapping = aes(features, 
                                              tsRNAs.df.subset$negLog10, 
                                              color=tsRNAs.df.subset$negLog10)) +
    geom_point() +
    ggtitle("DE Analysis - tsRNAs") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7)) +
    scale_color_gradient(low="blue", high="red") +
    #scale_y_continuous(trans='log2') +   # Change y axis to log scale
    scale_x_discrete(limits = (levels(tsRNAs.df.subset$negLog10))) +
    labs(colour = "-Log10 of padj", 
         x = "ncRNA/gene", 
         y = "-Log10 of padj", 
         subtitle = "Max number of features shown is 20"))
  dev.off()
  
  ### Plot genes and sno/miRNAs:
  pdf.width <- nrow(genes.df.subset)*0.2 + 3
  pdf(file = paste0(path.to.files, args[2], "_snomiRNAs-and-genes.high-DE-negLog10_padj.pdf"), width = pdf.width, height = 5)
  print(ggplot(data = genes.df.subset, mapping = aes(features, 
                                                genes.df.subset$negLog10, 
                                                color=genes.df.subset$negLog10)) +
    geom_point() +
    ggtitle("DE Analysis - Genes and sno/miRNAs") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7)) +
    scale_color_gradient(low="blue", high="red") +
    #scale_y_continuous(trans='log2') +   # Change y axis to log scale
    scale_x_discrete(limits = (levels(genes.df.subset$negLog10))) +
    labs(colour = "-Log10 of padj", 
         x = "ncRNA/gene", 
         y = "-Log10 of padj", 
         subtitle = "Max number of features shown is 20"))
  dev.off()

  }


#--------------------#
#    Run analysis    #
#--------------------#

DESeq2.function(myPath)

