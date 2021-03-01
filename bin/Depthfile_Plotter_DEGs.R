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
mean.cutoff <- as.integer(args[5])

### Print everything:
flag <- args[6]

### Check if the correct number of command line arguments were provided. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
}

### Load ggplot2
library(ggplot2)

### Read in files and get file basename
input1 <- read.table(args[1])
input2 <- read.table(args[2])
genelist <- read.table(args[3]) # File with list of genes to plot
genelist$V1 <- factor(genelist$V1, levels=unique(genelist$V1)) # Reassign the factor to stop reordering of factors when converting DF column to list
genelist <- genelist$V1 # Convert DF column to list

if (file.size(args[3]) == 0) { # Check if features file contains any lines
  message <- paste0("The file ", as.character(args[3]), " was empty. There were no features to plot")
  write.table(x = message, file = paste0(args[3], ".README"))
} else {
  #input3 <- read.table(args[3], sep = "\t")
  if (length(args)==7) {
    GTF <- read.table(args[7], sep = "\t")
  } 
  
  ### Get names of conditions from arg3 name for writing in plots
  file3 <- basename(args[4])
  conditions.for.plots <- strsplit(x = file3, "_Features")[[1]][1] #Must do double index to access the resulting list from strsplit
  condition1 <- strsplit(x = conditions.for.plots, "_vs_")[[1]][1]
  condition2 <- strsplit(x = conditions.for.plots, "_vs_")[[1]][2]
  
  ### Get features from files
  features <- input1$V1   # Group features by name
  featuresUnion <- union(features, features) # Get unique set of features
  
  ### Open a PDF for writing
  plot_list = list()
  
  ##### If all features are not to be plotted:
  if(flag == "no") {   
    counter = 0
    for(feature in genelist){
      if (any(grepl(feature,input1$V1))) { # Only generate plots for features in provided list
        # Do nothing and continue to generate plot etc.
      } else {
        next  # Skip this iteration and do not generate plot
      }
      if(counter == 20){ # Only print top 20 features
        break # exit for loop
      } else {
        counter = counter + 1
        if(grepl("?", feature)) {
          feature <- gsub("\\?", "\\\\?", feature)   # Allows script to grep for tRNA with ??? in name
        }
        ##### Format and subset dataframes for files 1 and 2
        ### file 1
        feature.full <- paste0("^", feature, "$", sep="")
        subset1 <- input1[grep(feature.full, input1$V1),] # \\ + \\b denote 'boundary' of word we are searching for. Stops partial matches. Equivalent to ^feature$
        subset1$coordinates <- 1:nrow(subset1)   # Make a new column with coordinates starting from 1
        subset1$Conditions <- condition1 # Create column with condition 1 in every row
        subset1.mean <- mean(subset1[,ncol(subset1)-3]) # Get the mean column at the end of the dataframe
        subset1.meancol.num <- ncol(subset1)-3  # Get the column index for mean of condition 1
        subset1.stdcol.num <- subset1.meancol.num + 1 # Get the column index for the std dev (mean column index + 1) 
        new.subset1 <- subset1[,c(1, subset1.meancol.num, subset1.stdcol.num, subset1.stdcol.num+1, subset1.stdcol.num+2)] # Subset dataframe to get specific rows
        names(new.subset1) <- c("feature", 
                                "mean",
                                "std",
                                "coord",
                                "Conditions") # Create column header names

        ### file 2 (see file 1 comments above for explanation)
        subset2 <- input2[grep(feature.full, input2$V1),]
        subset2$coordinates <- 1:nrow(subset2)
        subset2$Conditions <- condition2  
        subset2.mean <- mean(subset2[,ncol(subset2)-3])
        subset2.meancol.num <- ncol(subset2)-3
        subset2.stdcol.num <- subset2.meancol.num + 1
        new.subset2 <- subset2[,c(1, subset2.meancol.num, subset2.stdcol.num, subset2.stdcol.num+1, subset2.stdcol.num+2)]
        names(new.subset2) <- c("feature", 
                                "mean",
                                "std",
                                "coord",
                                "Conditions")
        
        ### get ncRNA gene names rather than IDs
        if (length(args)==7) {
          featureRows <- GTF[grep(feature, GTF$V9),]
          featureRows <- featureRows[1,]
          geneName <- as.character(sub(".*gene_name *(.*?) *; .*", "\\1", featureRows$V9))
          geneName <- paste0(geneName, " (", feature, ")")
        } else if (nchar(feature) == 3) {
          geneName <- paste0(feature, " (plot of multi-mapping reads where isoacceptor could not be identified)", sep="")
        } else {
          geneName <- feature
        }
        ### combine dataframes
        new.df2 <- rbind(new.subset1, new.subset2)

        ### Plot only if the mean is above the min (arg4)
        if(subset1.mean >= mean.cutoff | subset2.mean >= mean.cutoff){
          p = ggplot(new.df2, aes(x=coord, y=mean, fill=Conditions), show.legend = FALSE) +
            geom_ribbon(aes(ymin=mean-std, ymax=mean+std, colour=Conditions),
                        alpha=0.3) +
            geom_line(aes(linetype=Conditions)) + # Different linetypes
            #geom_line(size=0.8) +
            xlab('Nucleotide position') +
            ylab('Coverage (TPM)') +
            scale_fill_manual(values=c("red", "blue")) +
            scale_color_manual(values=c("red", "blue")) +
            labs(fill="Conditions", title = geneName)
          plot_list[[feature]] = p
      
          }
    
        }
    }
    pdf(args[4])
    for(feature in genelist) {
      capture.output(print(plot_list[[feature]])) # Capture.output suppresses print output 'NULL'
    }
    invisible(dev.off()) # Suppress message 'NULL'
    
  } else {
    for(feature in featuresUnion) {
      if(grepl("?", feature)) {
        feature <- gsub("\\?", "\\\\?", feature)   # Allows script to grep for tRNA with ??? in name
      }
      
      ##### Format and subset dataframes for files 1 and 2
      ### file 1
      feature.full <- paste0("^", feature, "$", sep="")
      subset1 <- input1[grep(feature.full, input1$V1),] # \\ + \\b denote 'boundary' of word we are searching for. Stops partial matches. Equivalent to ^feature$
      subset1$coordinates <- 1:nrow(subset1)   # Make a new column with coordinates starting from 1
      subset1$Conditions <- condition1 # Create column with condition 1 in every row
      subset1.mean <- mean(subset1[,ncol(subset1)-3]) # Get the mean column at the end of the dataframe
      subset1.meancol.num <- ncol(subset1)-3  # Get the column index for mean of condition 1
      subset1.stdcol.num <- subset1.meancol.num + 1 # Get the column index for the std dev (mean column index + 1) 
      new.subset1 <- subset1[,c(1, subset1.meancol.num, subset1.stdcol.num, subset1.stdcol.num+1, subset1.stdcol.num+2)] # Subset dataframe to get specific rows
      names(new.subset1) <- c("feature", 
                              "mean",
                              "std",
                              "coord",
                              "Conditions") # Create column header names
      
      ### file 2 (see file 1 comments above for explanation)
      subset2 <- input2[grep(feature.full, input2$V1),]
      subset2$coordinates <- 1:nrow(subset2)
      subset2$Conditions <- condition2  
      subset2.mean <- mean(subset2[,ncol(subset2)-3])
      subset2.meancol.num <- ncol(subset2)-3
      subset2.stdcol.num <- subset2.meancol.num + 1
      new.subset2 <- subset2[,c(1, subset2.meancol.num, subset2.stdcol.num, subset2.stdcol.num+1, subset2.stdcol.num+2)]
      names(new.subset2) <- c("feature", 
                              "mean",
                              "std",
                              "coord",
                              "Conditions")
      
      ### get ncRNA gene names rather than IDs
      if (length(args)==7) {
        featureRows <- GTF[grep(feature, GTF$V9),]
        featureRows <- featureRows[1,]
        geneName <- as.character(sub(".*gene_name *(.*?) *; .*", "\\1", featureRows$V9))
        geneName <- paste0(geneName, " (", feature, ")")
      } else if (nchar(feature) == 3) {
        geneName <- paste0(feature, " (plot of multi-mapping reads where isoacceptor could not be identified)", sep="")
      } else {
        geneName <- feature
      }
      ### combine dataframes
      new.df2 <- rbind(new.subset1, new.subset2)
      
      ### Plot only if the mean is above the min (arg4)
      if(subset1.mean >= mean.cutoff | subset2.mean >= mean.cutoff){
        p = ggplot(new.df2, aes(x=coord, y=mean, fill=Conditions), show.legend = FALSE) +
          geom_ribbon(aes(ymin=mean-std, ymax=mean+std, colour=Conditions),
                      alpha=0.3) +
          geom_line(aes(linetype=Conditions)) + # Different linetypes
          #geom_line(size=0.8) +
          xlab('Nucleotide position') +
          ylab('Coverage (TPM)') +
          scale_fill_manual(values=c("red", "blue")) +
          scale_color_manual(values=c("red", "blue")) +
          labs(fill="Conditions", title = geneName)
        plot_list[[feature]] = p
      }
  }
  pdf(args[4])
  for(feature in featuresUnion) {
    capture.output(print(plot_list[[feature]])) # Capture.output suppresses print output 'NULL'
  }
  invisible(dev.off()) # Suppress message 'NULL'
  }
}

