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
mean.cutoff <- args[4]

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
}

### Load ggplot2
library(ggplot2)

### Read in files and get file basename
input1 <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/NewOutput2/Results/Merged_mean-std1.txt.csv")
file1 <- basename("/home/paul/Documents/Pipelines/tirna-pipeline/NewOutput2/Results/Merged_mean-std1.txt.csv")
input2 <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/NewOutput2/Results/Merged_mean-std2.txt.csv")
file2 <- basename("/home/paul/Documents/Pipelines/tirna-pipeline/NewOutput2/Results/Merged_mean-std2.txt.csv")

### Get features from files
features <- input1$V1   # Group features by name
featuresUnion <- union(features, features) # Get unique set of features

### Open a PDF for writing
plot_list = list()

#pdf(args[3])
for(feature in featuresUnion) {
  
  ### file 1
  subset1 <- input1[grep(feature, input1$V1),]
  subset1.NoFlanks <- tail(subset1, -10)
  subset1.NoFlanks <- head(subset1.NoFlanks, -10)
  subset1.NoFlanks$coordinates <- 1:nrow(subset1.NoFlanks)   # Make a new column with coordinates starting from 1
  subset1.NoFlanks$Conditions <- file1
  subset1.mean <- mean(subset1.NoFlanks$V10)
  
  ### file 2
  subset2 <- input2[grep(feature, input2$V1),]
  subset2.NoFlanks <- tail(subset2, -10)
  subset2.NoFlanks <- head(subset2.NoFlanks, -10)
  subset2.NoFlanks$coordinates <- 1:nrow(subset2.NoFlanks)   # Make a new column with coordinates starting from 1
  subset2.NoFlanks$Conditions <- file2   
  subset2.mean <- mean(subset2.NoFlanks$V10)
  
  
  ### combine dataframes
  new.df <- rbind(subset1.NoFlanks, subset2.NoFlanks)
  
  ### Plot only if the mean is above the min (set at the beginning of script)
  if(subset.mean>mean.cutoff){

  p = ggplot(new.df, aes(x=coordinates, y=V10, fill=Conditions), show.legend = FALSE) +
      geom_ribbon(aes(ymin=V10-V11, ymax=V10+V11, colour=Conditions),
                  alpha=0.3) +
      geom_line(size=0.8) +
      xlab('Nucleotide position') +
      ylab('Read coverage') +
      scale_fill_manual(values=c("red", "blue")) +
      scale_color_manual(values=c("red", "blue")) +
      labs(fill="Conditions", title = feature)
  plot_list[[feature]] = p
  }
}

#pdf("/home/paul/Documents/Pipelines/tirna-pipeline/NewOutput2/Results/Output.pdf")
pdf(args[3])
for(feature in featuresUnion) {
  print(plot_list[[feature]])
}
dev.off()


