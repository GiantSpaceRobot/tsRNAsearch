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
mean.cutoff <- args[4]

### Check if the correct number of command line arguments were provided. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
}

### Load ggplot2
library(ggplot2)

### Read in files and get file basename
#input1 <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/NewOutput2/Results/Data/condition1_concatenated_mean_stdev.genomecov")
#input2 <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/NewOutput2/Results/Data/condition2_concatenated_mean_stdev.genomecov")
#mean.cutoff = 20
input1 <- read.table(args[1])
input2 <- read.table(args[2])

### Get names of conditions from arg3 name for writing in plots
file3 <- basename(args[3])
#file3 <- "condition1_vs_condition2_DEGs.pdf"
conditions.for.plots <- strsplit(x = file3, "_DEGs")[[1]][1] #Must do double index to access the resulting list from strsplit
condition1 <- strsplit(x = conditions.for.plots, "_vs_")[[1]][1]
condition2 <- strsplit(x = conditions.for.plots, "_vs_")[[1]][2]

### Get features from files
features <- input1$V1   # Group features by name
featuresUnion <- union(features, features) # Get unique set of features

### Open a PDF for writing
plot_list = list()

#pdf("/home/paul/Documents/Pipelines/tirna-pipeline/Output_NoStdev.pdf")

pdf(args[3])
for(feature in featuresUnion) {
  
  ### file 1
  subset1 <- input1[grep(feature, input1$V1),]
  subset1.NoFlanks <- tail(subset1, -10)
  subset1.NoFlanks <- head(subset1.NoFlanks, -10)
  subset1.NoFlanks$coordinates <- 1:nrow(subset1.NoFlanks)   # Make a new column with coordinates starting from 1
  subset1.NoFlanks$Conditions <- condition1
  subset1.mean <- mean(subset1.NoFlanks$V10)
  
  ### file 2
  subset2 <- input2[grep(feature, input2$V1),]
  subset2.NoFlanks <- tail(subset2, -10)
  subset2.NoFlanks <- head(subset2.NoFlanks, -10)
  subset2.NoFlanks$coordinates <- 1:nrow(subset2.NoFlanks)   # Make a new column with coordinates starting from 1
  subset2.NoFlanks$Conditions <- condition2  
  subset2.mean <- mean(subset2.NoFlanks$V10)
  
  ### combine dataframes
  new.df <- rbind(subset1.NoFlanks, subset2.NoFlanks)
  new.df2 <- cbind(subset1.NoFlanks, subset2.NoFlanks)
  cols <- seq(1, length(new.df2))
  cols <- paste0("V", cols)
  colnames(new.df2) <- c(cols)
  
  ### Plot only if the mean is above the min (arg4)
  if(subset1.mean > mean.cutoff | subset2.mean > mean.cutoff){
  #p = ggplot(new.df2) + 
  #        geom_line(aes(x=V12, y=V10),
  #            color = "red") +
  #        geom_line(aes(x=V12, y=V23),
  #          color = "blue") +
  #        xlab('Nucleotide position') +
  #        ylab('Read coverage') +
  #        labs(fill="Conditions", title = feature) +
  #        annotate("text", x = 2, y=-100, label = condition1, color = "red") + 
  #        annotate("text", x = 20, y=-100, label = condition2, color = "blue") 
  min.y <- min(c(new.df2$V10, new.df2$V23))
  max.y <- max(c(new.df2$V10, new.df2$V23))
  plot(x = new.df2$V12, y = new.df2$V10, 
       ylim = c(min.y, max.y),
       #log = "y", 
       type = "l", 
       col = "red", 
       main=feature, 
       xlab="Nucleotide position", 
       ylab="Read coverage") 
  #points(x = new.df2$V12, y = new.df2$V10 , type = "l", col = 1) +
  points(x = new.df2$V12, y = new.df2$V23, 
         #log = "y", 
         type = "l", 
         col = "blue") 
  #legend("bottomleft", legend = (c(condition1, condition2)), border = "black", bg = "white", col=c("red","blue"), pch=1)
  legend("bottomright", legend = (c(condition1, condition2)), border = "black", bg = "white", col=c("red","blue"), pch=1)
  #plot_list[[feature]] = p
  }
}
dev.off()


#pdf("/home/paul/Documents/Pipelines/tirna-pipeline/Output_NoStdev.pdf")
#pdf(args[3])
#for(feature in featuresUnion) {
#  print(plot_list[[feature]])
#}
#dev.off()


