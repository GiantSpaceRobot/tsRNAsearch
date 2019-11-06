#!/usr/bin/env Rscript

## Script name: DESeq2_tsRNAsearch.R
##
## Purpose of script: Carry out differential gene expression analysis of ncRNAs from small RNA-seq data
##
## Author: Dr Paul Donovan
##
## Date Created: 2-1-2019 (2nd of Jan 2019)
##
## Copyright: MIT license
##
## Email: pauldonovan@rcsi.com

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  message("No arguments provided. Defaulting to condition 1 vs condition 2 assuming three replicates each.")
} else if (length(args)>4) {
  stop("Error: Too many command line arguments. Quitting.")
} else if (length(args)==1) {
  stop("Error: Require three command line arguments at least: path-to-files, output-prefix, ncRNA GTF, [layout-file.csv]. Quitting.")
} else if (length(args)==2) {
  stop("Error: Require three command line arguments at least: path-to-files, output-prefix, ncRNA GTF, [layout-file.csv]. Quitting.")
} else if (length(args)==3) {             # If only a filepath was provided:
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
} else if (length(args)==4) {     # If a .csv and file path were provided:
  experiment <- args[4]
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
  file.CSV <- read.csv(args[4], header=FALSE)
  #print(file.CSV)
  lvls.df <- as.data.frame(table(file.CSV$V2))
  #file.CSV <- read.csv("/home/paul/Documents/Data/EpimiRNA_Data/Comparisons/Mouse_IAKA_KA2w-vs-KA2w-ctrl_Layout.csv", header = F)
  #table(file.CSV)
  #print(lvls.df)
  #arg2 <- "KA_2w_vs_control"
  my.conditions <- unlist(strsplit(args[2], "_vs_"))
  #my.conditions <- unlist(strsplit(arg2, "_vs_"))
  Condition1 <- my.conditions[1]
  Condition2 <- my.conditions[2]
  ReplicateNumber1 <- lvls.df[grep(Condition1, lvls.df$Var1),][1,2]
  ReplicateNumber2 <- lvls.df[grep(Condition2, lvls.df$Var1),][1,2]
  #condition2reps <- paste0(Condition2, 1:ReplicateNumber2)
  #class(ReplicateNumber1)
  #ReplicateNumber1 <- lvls.df[1,2]
  #Condition1 <- toString(lvls.df[1,1])
  #ReplicateNumber2 <- lvls.df[2,2]
  #Condition2 <- toString(lvls.df[2,1])
  #print(Condition1)
  #print(Condition2)
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
ResultsFile <- paste0(Condition1,"_vs_",Condition2)

### Directory for writing DESeq2 results to file
if (file.exists(paste0(myPath, "DE_Results/DESeq2"))){
} else {
  dir.create(paste0(myPath, "DE_Results/DESeq2"))
}

#-----------------#
#    Libraries    #
#-----------------#

# If libraries not installed, install them
#source("http://bioconductor.org/biocLite.R")
#biocLite(DESeq2)
#library(devtools)
#devtools::install_github('kevinblighe/EnhancedVolcano')

suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(reshape2))

### Reading in gene mapping file
GTF <- read.table(args[3], sep = "\t") # Read in GTF for name conversion
GTF$NewNames <- paste0(GTF$V1, " (", as.character(sub(".*gene_name *(.*?) *; .*", "\\1", GTF$V9)), ")") # Add new column with gene names
GTF <- GTF[!duplicated(GTF$NewNames),]  # Remove duplicates based on NewNames

#------------------------------#
#       Define functions       #
#------------------------------#

### Function to carry out DESeq2 analysis on provided files
DESeq2.function <- function(path.to.files){
  
  ### Create checkpoints (for simple degbugging) 
  print("Checkpoint 1")

  ### Read files
  path.to.files <- myPath
  file.names <- dir(path.to.files, pattern =".count")
  #print(file.names)
  cDataAll <- NULL
  for (i in 1:length(file.names)){
    full.path <- paste0(path.to.files, file.names[i])
    file <- read.table(full.path, header = TRUE)
    cDataAll <- cbind(cDataAll, file[,2])
  }
  #Condition1 <- "CytC"
  #Condition2 <- "CytC-ctrl"
  condition1reps <- paste0(Condition1, 1:ReplicateNumber1)
  condition2reps <- paste0(Condition2, 1:ReplicateNumber2)
  newcolnames <- c(condition1reps, condition2reps) # Create column names using the provided condition names
  #print(head(cDataAll))
  colnames(cDataAll) <- (newcolnames)
  #print(head(cDataAll))
  #colnames(cDataAll) <- (file.names)
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
  if((ReplicateNumber1==1) & (ReplicateNumber2==1)) {    # If replicate number = 1, stop analysis
    print("    Single replicate analysis, skipping formal DESeq2 analysis and carrying out single Log2 fold comparison instead...")
    log2FC.cutoff <- 1.5  # Log2FC cut-off can be changed here
    log2.df <- log2(tpm)  # Log2 transformation
    log2.df[which(!is.finite(log2.df))] <- 0 # Convert all Inf/-Inf/NA to 0
    log2FC <- (log2.df[,1] - log2.df[,2]) # Get Log2 fold change
    log2FC.df <- data.frame(log2.df[,1], log2.df[,2], log2FC)
    up <- (log2FC.df[!is.na(log2FC.df$log2FC) & log2FC.df$log2FC >= log2FC.cutoff, ])    
    down <- (log2FC.df[!is.na(log2FC.df$log2FC) & log2FC.df$log2FC <= -log2FC.cutoff, ]) 
    colnames(log2FC.df) <- c(Condition1, Condition2, "Log2FC")
    colnames(up) <- c(Condition1, Condition2, "Log2FC")
    colnames(down) <- c(Condition1, Condition2, "Log2FC")
    write.csv(log2FC.df, file=paste0(path.to.files, "DE_Results/DESeq2/", ResultsFile, "_DESeq2-output.csv"))
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
      
  } else {   # If replicate numbers are greater than 1:
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
  sampleDists <- data.frame(sampleDists)
  colnames(sampleDists) <- gsub(x = colnames(sampleDists), pattern = ".collapsed.all.features.count", replacement = "") # Remove string from colnames
  rownames(sampleDists) <- gsub(x = rownames(sampleDists), pattern = ".collapsed.all.features.count", replacement = "") # Remove string from rownames 
  #(mycols <- brewer.pal(length(file.names), "Paired")[1:length(unique(file.names))])
  pdf(paste0(path.to.files, "DE_Results/", ResultsFile, "_Distance-Matrix.pdf"),
      width=12,height=12)
  par(mar=c(6,4,4,5)+0.1) 
  heatmap.2(as.matrix(sampleDists), key=F, trace="none",
            col=colorpanel(100, "black", "white"),
            #ColSideColors=mycols[file.names], 
            #RowSideColors=mycols[file.names],
            #cexRow = 0.8,
            #cexCol = 0.8,
            margins=c(12,10),
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
  names(x) <- gsub(x = names(x), pattern = ".collapsed.all.features.count", replacement = "") # Remove string from names
  #names(y) <- gsub(x = names(y), pattern = ".collapsed.all.features.count", replacement = "") # Remove string from names
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
         inset=c(-0.2,0), # -0.5,0 to push labels to the right 
         pch=20,
         col=colLabel,
         cex = 0.7, 
         legend=names(x),
         xpd = TRUE,
         bty = "n",
         lty=NULL)
  garbage <- dev.off()
  #par(mar=c(5, 4, 4, 2) + 0.1)
  
  ### checkpoint
  print("Checkpoint 8")

  ### Convert ncRNA names into gene names
  print("Converting gene IDs to gene names")
  results.DF <- as.data.frame(res)
  results.DF$features <- rownames(results.DF) # Add column with gene IDs
  # Map names in GTF to gene ENSG IDs
  for( row in seq_len( nrow(GTF))) {  
    results.DF$mapped[ substr( results.DF$features, 0, nchar(as.character(GTF$V1[row]))) == GTF$V1[row] ] <- GTF$NewNames[row]
  }
  results.DF$mapped <- as.character(results.DF$mapped) # Convert from factor to character
  results.DF$features <- as.character(results.DF$features) # Convert from factor to character
  results.DF$final <- ifelse(is.na(results.DF$mapped), results.DF$features, results.DF$mapped) # Create new column using new name where possible
  resultsDF <- subset(results.DF, select = -c(features, mapped)) # Remove unneeded columns
  rownames(resultsDF) <- resultsDF$final # Replace rownames with new names
  resultsDF <- subset(resultsDF, select = -c(final)) # Remove unneeded columns
  
  ### checkpoint
  print("Checkpoint 9")
  
  write.csv(resultsDF, file=paste0(path.to.files, "DE_Results/DESeq2/", ResultsFile, "_DESeq2-output.csv"))
  
  up <- (resultsDF[!is.na(resultsDF$padj) & resultsDF$padj <= 0.05 &
               resultsDF$log2FoldChange >= 0.5, ])
  down <- (resultsDF[!is.na(resultsDF$padj) & resultsDF$padj <= 0.05 &    
                 resultsDF$log2FoldChange <= -0.5, ]) 
  
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
  print("Checkpoint 10")
  
  ### Generate data for plots
  newdata.subset <- rbind(up, down)
  newdata.subset$features <- rownames(newdata.subset)
  newdata.subset$negLog10 <- -log10(newdata.subset$padj)
  ncRNAs.df <- subset(newdata.subset, startsWith(as.character(features), "ENS"))
  tsRNAs.df <- newdata.subset[ !(newdata.subset$features %in% ncRNAs.df$features), ] # newdata.susbet minus all rows in ncRNAs.df (to get tsRNAs.df)
  #tsRNAs.df <- subset(newdata.subset, startsWith(as.character(features), "chr"))
  # If there are more than 20 features, show top 20
  if(nrow(tsRNAs.df) > 20){
    tsRNAs.df.subset <- head(tsRNAs.df, n = 20)
  } else {
    tsRNAs.df.subset <- tsRNAs.df
  }

  # If there are more than 20 features, show top 20
  if(nrow(ncRNAs.df) > 20){
    ncRNAs.df.subset <- head(ncRNAs.df, n = 20)
  } else {
    ncRNAs.df.subset <- ncRNAs.df
  }
  ### replace Inf values (extremely low adj p-value) with -Log10 of 300
  tsRNAs.df.subset <- data.frame(lapply(tsRNAs.df.subset, function(x) {gsub(Inf, "300", x)}))
  tsRNAs.df.subset$negLog10 <- as.numeric(levels(tsRNAs.df.subset$negLog10))[tsRNAs.df.subset$negLog10] #Convert factor type to numeric
  tsRNAs.df.subset$features <- factor(tsRNAs.df.subset$features, levels = tsRNAs.df.subset$features[order(tsRNAs.df.subset$negLog10)])
  ncRNAs.df.subset <- data.frame(lapply(ncRNAs.df.subset, function(x) {gsub(Inf, "300", x)}))
  ncRNAs.df.subset$negLog10 <- as.numeric(levels(ncRNAs.df.subset$negLog10))[ncRNAs.df.subset$negLog10] #Convert factor type to numeric
  ncRNAs.df.subset$features <- factor(ncRNAs.df.subset$features, levels = ncRNAs.df.subset$features[order(ncRNAs.df.subset$negLog10)])
  
  ### checkpoint
  print("Checkpoint 11")
  
  pdf.width <- nrow(tsRNAs.df.subset)*0.2 + 3
  pdf(file = paste0(path.to.files, args[2], "_tsRNAs.high-DE-negLog10_padj.pdf"), width = pdf.width, height = 5)
  print(ggplot(data = tsRNAs.df.subset, mapping = aes(features, 
                                              tsRNAs.df.subset$negLog10, 
                                              color=tsRNAs.df.subset$negLog10)) +
    geom_point() +
    ggtitle("DE Analysis - tRNAs") +
    theme(axis.text.x = element_text(size=7)) + # Move x axis label down
    scale_color_gradient(low="blue", high="red") +
    #scale_y_continuous(trans='log2') +   # Change y axis to log scale
    scale_x_discrete(limits = (levels(tsRNAs.df.subset$negLog10))) +
    coord_flip() +
    labs(colour = "-Log10\n   of \n  padj", 
         x = "tRNAs", 
         y = "-Log10 of padj", 
         subtitle = "Max number of features shown is 20"))
  dev.off()
  
  ### Plot ncRNAs:
  pdf.width <- nrow(ncRNAs.df.subset)*0.2 + 3
  pdf(file = paste0(path.to.files, args[2], "_ncRNAs.high-DE-negLog10_padj.pdf"), width = pdf.width, height = 5)
  print(ggplot(data = ncRNAs.df.subset, mapping = aes(features, 
                                                ncRNAs.df.subset$negLog10, 
                                                color=ncRNAs.df.subset$negLog10)) +
    geom_point() +
    ggtitle("DE Analysis - ncRNAs") +
    theme(axis.text.x = element_text(size=7)) +
    scale_color_gradient(low="blue", high="red") +
    #scale_y_continuous(trans='log2') +   # Change y axis to log scale
    scale_x_discrete(limits = (levels(ncRNAs.df.subset$negLog10))) +
    coord_flip() +
    labs(colour = "-Log10\n   of \n  padj", 
         x = "ncRNAs", 
         y = "-Log10 of padj", 
         subtitle = "Max number of features shown is 20"))
  dev.off()
  
  ### checkpoint
  print("Checkpoint 12")
  
  ### Create Volcano plot using EnhancedVolcano fro Kevin Blighe
  volcano <- EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  xlim = c(-5, 8),  
                  pCutoff = 0.05, 
                  FCcutoff = 0.5,
                  labSize = 3.0)
  volcano + 
    labs(subtitle = "") # Remove automatic subtitle
  ggsave(filename = paste0(path.to.files, args[2], "_VolcanoPlot.pdf"), plot = volcano) # Save plot using ggplot2 ggsave (error occured using normal R PDF save)
  
  ### checkpoint
  print("Checkpoint 13")
  
  ### Create barplots for all features identified
  my.levels <- as.character(levels.default(groups)) # Get names of levels
  #print(my.levels)
  #my.levels <- c(Condition1, Condition2)
  level1 <- my.levels[1] # get name of level 1
  level2 <- my.levels[2] # get name of level 2
  #print(head(cDataAll))
  #library(dplyr)
  #level1 <- "CytC"
  #level2 <- "CytC-ctrl"
  #cDataAll <- data.frame("CytC1"=0,
  #                 "CytC2"=0, 
  #                 "CytC-ctrl1"=0,
  #                 "CytC-ctrl2"=0) 
  
  #level1.df <- select(data.frame(cDataAll), contains(level1)) # Subset main count matrix using condition 1
  level1.df <- data.frame(cDataAll) %>% select(c(1:ReplicateNumber1)) # Subset main count matrix to get condition 1
  #level2.df <- select(data.frame(cDataAll), contains(level2)) # Subset main count matrix using condition 2
  level2.df <- data.frame(cDataAll) %>% select(c((ReplicateNumber1+1):(ReplicateNumber1+ReplicateNumber2))) # Subset main count matrix to get condition 2
  #print("first")
  #print(head(level1.df))
  #print("second")
  #print(head(level2.df))
  level1.df$total.raw <- rowSums(level1.df) #Calculate total read count for condition 1
  level2.df$total.raw <- rowSums(level2.df) #Calculate total read count for condition 1
  level1.df$total.rpm <- level1.df$total.raw/(sum(level1.df$total.raw)/1000000) # Calculate reads per million
  level2.df$total.rpm <- level2.df$total.raw/(sum(level2.df$total.raw)/1000000) # Calculate reads per million
  rpm.plot <- create.barplot(level1.df$total.rpm, level2.df$total.rpm, "RPM", rownames(level1.df), level1, level2) # Call function
  raw.plot <- create.barplot(level1.df$total.raw, level2.df$total.raw, "Raw Read Count", rownames(level1.df), level1, level2) # Call function
  ggsave(filename = paste0(path.to.files, args[2], "_BarPlot_RPM-normalised.pdf"), plot = rpm.plot) # Save plot using ggplot2 ggsave 
  ggsave(filename = paste0(path.to.files, args[2], "_BarPlot_Raw-readcounts.pdf"), plot = raw.plot) # Save plot using ggplot2 ggsave 
  
  } # Checkpoint 3 - single replicate analysis if statement
} # End of function

### Function to create barplots
create.barplot <- function(column.cond1, column.cond2, normalisation.method, my.rownames, name1, name2) {
  subtitle <- ifelse(normalisation.method == "RPM", "Normalised to reads per million (RPM)", "Raw read counts displayed") # Create subtitle for plot
  #print(column.cond1)
  #print(column.cond2)
  barplot.df <- data.frame(cbind(column.cond1, column.cond2))
  names(barplot.df)[1] <- "cond1" #level1 # Add correct colname to DF
  names(barplot.df)[2] <- "cond2" #level2 # Add correct colname to DF
  rownames(barplot.df) <- my.rownames
  barplot.df$features <- my.rownames
  GTF$ncRNA.type <- gsub(".*\\gene_biotype ", "", GTF$V9) # Make a new column with ncRNA type (remove everything before gene_biotype in col 9)
  GTF$ncRNA.type <- gsub(";.*", "", GTF$ncRNA.type) # Replace the semicolon at the end of string
  for(row in seq_len(nrow(GTF))) {  
    barplot.df$mapped[ substr( barplot.df$features, 0, nchar(as.character(GTF$V1[row]))) == GTF$V1[row] ] <- GTF$ncRNA.type[row]
  }
  barplot.df$mapped <- as.character(barplot.df$mapped) # Convert from factor to character
  #print(head(barplot.df))
  ### Assign ncRNA labels for as many ncRNA types as possible
  for(row in seq_len( nrow(barplot.df))) {
    if(is.na(barplot.df$mapped[row])) {   # If no ncRNA type has been assigned yet:
      if(startsWith(barplot.df$features[row], "ENS") == TRUE) {  # If the feature has an 'ENS' ID, it is a gene or unknown feature. Ignore for now.
        barplot.df$mapped[row] <- "other (ncRNA)"
      } else { # If the feature has no ncRNA type yet and does not start with 'ENS', it is a tRNA
        barplot.df$mapped[row] <- "tRNA"
      }
    } 
  } ### End of ncRNA label assignment loop
  #barplot.df <- barplot.df[complete.cases(barplot.df),]
  barplot.df$mapped <- as.factor(barplot.df$mapped)
  grouped.df <- barplot.df %>% group_by(mapped) %>% summarise(cond1.sum = sum(cond1), cond2.sum = sum(cond2))
  names(grouped.df)[1] <- "Feature" # Add correct colname to DF
  names(grouped.df)[2] <- name1 # Add correct colname to DF
  names(grouped.df)[3] <- name2 # Add correct colname to DF
  #print(head(grouped.df))
  melted.df <- melt(data = grouped.df, id.vars = "Feature")
  df.cols <- length(unique(melted.df$Feature)) # Expand color palette
  mycolors <- colorRampPalette(brewer.pal(11, "RdYlBu"))(df.cols) # Expand RdYlBu to # of colours needed based on DF
  results.barplot <- ggplot(data = melted.df, mapping = aes(variable, value, fill = Feature)) +
    geom_bar(stat = "identity", width = 0.5) + 
    ggtitle(paste0("Features identified in ", name1, " vs ", name2)) +
    scale_fill_manual(values = mycolors) +
    #theme_minimal() +
    #coord_flip() +
    theme(legend.position = "bottom") +
    labs(x = "", 
         y = normalisation.method, 
         subtitle = subtitle)
  #results.barplot <- ggplot(data = melted.df, mapping = aes(variable, value, fill = Feature)) +
  #  geom_bar(stat = "identity", width = 0.5) + 
  #  ggtitle(paste0("Features identified in ", name1, " vs ", name2)) +
  #  scale_fill_brewer(palette = "RdYlBu") +
  #  #theme_minimal() +
  #  #coord_flip() +
  #  labs(x = "", 
  #       y = normalisation.method, 
  #       subtitle = subtitle)
  return(results.barplot)
}

#--------------------#
#    Run analysis    #
#--------------------#

DESeq2.function(myPath)

