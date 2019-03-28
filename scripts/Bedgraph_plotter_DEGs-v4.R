#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script makes coverage plots
### 
###-------------------------------------------------------------------------------

### Use command line arguments
args = commandArgs(trailingOnly=TRUE)

### Change this if you want only to plot features with a higher mean coverage
### (e.g. mean coverage of 100 reads is 'mean.cutoff <- 100')
#mean.cutoff <- 0
mean.cutoff <- as.integer(args[5])

### Check if the correct number of command line arguments were provided. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
}

### Load ggplot2
library(ggplot2)

### Read in files and get file basename
#input1 <- read.table("/home/paul/Documents/Pipelines/tsRNAsearch/subset-fake/Results/Data/Intermediate-files/condition1_concatenated_mean_stdev.tiRNA.depth")
#input2 <- read.table("/home/paul/Documents/Pipelines/tsRNAsearch/subset-fake/Results/Data/Intermediate-files/condition2_concatenated_mean_stdev.tiRNA.depth")
#input3 <- read.table("/home/paul/Documents/Pipelines/tsRNAsearch/subset-fake/Results/Data/DE_Results/DESeq2/DEGs_names-only.txt")
#mean.cutoff = 0

input1 <- read.table(args[1])
input2 <- read.table(args[2])
if (file.size(args[3]) == 0) { # Check if features file contains any lines
  message <- paste0("The file ", as.character(args[3]), " was empty. There were no features to plot")
  write.table(x = message, file = paste0(args[3], ".README"))
} else {

  input3 <- read.table(args[3])
  
  #if (nrow(input3) < 1) {
  #  message <- paste0("The file ", args[3], " was empty. There were no features to plot")
  #  write.table(x = message, file = paste0(args[3], ".txt"))
  #} else {
  
  if (length(args)==6) {
    GTF <- read.table(args[6], sep = "\t")
    #GTF <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/DBs/hg19-snomiRNA.gtf", sep = "\t")
  } 
  
  ### Get names of conditions from arg3 name for writing in plots
  #file3 <- "Angio_vs_Vehic_Features_DiffExprtiRNAs.pdf"
  file3 <- basename(args[4])
  conditions.for.plots <- strsplit(x = file3, "_Features")[[1]][1] #Must do double index to access the resulting list from strsplit
  condition1 <- strsplit(x = conditions.for.plots, "_vs_")[[1]][1]
  condition2 <- strsplit(x = conditions.for.plots, "_vs_")[[1]][2]
  
  ### Get features from files
  features <- input1$V1   # Group features by name
  featuresUnion <- union(features, features) # Get unique set of features
  
  ### Open a PDF for writing
  plot_list = list()
  
  #pdf(args[3])
  for(feature in featuresUnion) {
    if (any(grepl(feature,input3$V1))) { # Only generate plots for differentially expressed genes
      # Do nothing and continue to generate plot etc.
    } else next  # Skip this iteration and do not generate plot
    
    ### file 1
    #feature <- "ValCAC"
    subset1 <- input1[grep(feature, input1$V1),]
    #subset1.NoFlanks <- tail(subset1, -10)
    #subset1.NoFlanks <- head(subset1.NoFlanks, -10)
    subset1.NoFlanks <- subset1
    subset1.NoFlanks$coordinates <- 1:nrow(subset1.NoFlanks)   # Make a new column with coordinates starting from 1
    subset1.NoFlanks$Conditions <- condition1
    subset1.mean <- mean(subset1.NoFlanks[,ncol(subset1.NoFlanks)-3]) # Get the mean column at the end of the dataframe
    subset1.meancol.num <- ncol(subset1.NoFlanks)-3
    subset1.stdcol.num <- ncol(subset1.NoFlanks)-2
    new.subset1 <- subset1.NoFlanks[,c(1, subset1.meancol.num, subset1.stdcol.num, subset1.stdcol.num+1, subset1.stdcol.num+2)]
    names(new.subset1) <- c("feature", 
                            "mean",
                            "std",
                            "coord",
                            "Conditions")
    
    ### file 2
    subset2 <- input2[grep(feature, input2$V1),]
    #subset2.NoFlanks <- tail(subset2, -10)
    #subset2.NoFlanks <- head(subset2.NoFlanks, -10)
    subset2.NoFlanks <- subset2
    subset2.NoFlanks$coordinates <- 1:nrow(subset2.NoFlanks)   # Make a new column with coordinates starting from 1
    subset2.NoFlanks$Conditions <- condition2  
    subset2.mean <- mean(subset2.NoFlanks[,ncol(subset2.NoFlanks)-3])
    subset2.meancol.num <- ncol(subset2.NoFlanks)-3
    subset2.stdcol.num <- ncol(subset2.NoFlanks)-2
    new.subset2 <- subset2.NoFlanks[,c(1, subset2.meancol.num, subset2.stdcol.num, subset2.stdcol.num+1, subset2.stdcol.num+2)]
    names(new.subset2) <- c("feature", 
                            "mean",
                            "std",
                            "coord",
                            "Conditions")
    
    #geneName <- feature
    ### get sno/miRNA gene names rather than IDs
    if (length(args)==6) {
      featureRows <- GTF[grep(feature, GTF$V9),]
      featureRows <- featureRows[1,]
      geneName <- as.character(sub(".*gene_name *(.*?) *; gene_source.*", "\\1", featureRows$V9))
      geneName <- paste0(geneName, " (", feature, ")")
    } else {
      geneName <- feature
    }
    
    ### combine dataframes
    new.df2 <- rbind(new.subset1, new.subset2)
    
    ### Plot only if the mean is above the min (arg4)
    if(subset1.mean > mean.cutoff | subset2.mean > mean.cutoff){
    p = ggplot(new.df2, aes(x=coord, y=mean, fill=Conditions), show.legend = FALSE) +
        geom_ribbon(aes(ymin=mean-std, ymax=mean+std, colour=Conditions),
                    alpha=0.3) +
        geom_line(size=0.8) +
        xlab('Nucleotide position') +
        ylab('Coverage (reads per million)') +
        scale_fill_manual(values=c("red", "blue")) +
        scale_color_manual(values=c("red", "blue")) +
        labs(fill="Conditions", title = geneName)
    plot_list[[feature]] = p
    }
  }
  
  #pdf("/home/paul/Documents/Pipelines/tirna-pipeline/Output2.pdf")
  pdf(args[4])
  for(feature in featuresUnion) {
    print(plot_list[[feature]])
  }
  garbage <- dev.off() # Suppress message 'NULL'
}
