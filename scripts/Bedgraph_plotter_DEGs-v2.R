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
mean.cutoff <- as.integer(args[4])

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

if (length(args)==4) {
  GTF <- read.table(args[4], sep = "\t")
  #GTF <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/DBs/hg19-snomiRNA.gtf", sep = "\t")
} 

### Get names of conditions from arg3 name for writing in plots
file3 <- basename(args[3])
conditions.for.plots <- strsplit(x = file3, "_DEGs")[[1]][1] #Must do double index to access the resulting list from strsplit
condition1 <- strsplit(x = conditions.for.plots, "_vs_")[[1]][1]
condition2 <- strsplit(x = conditions.for.plots, "_vs_")[[1]][2]

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
  subset1.NoFlanks$Conditions <- condition1
  subset1.mean <- mean(subset1.NoFlanks$V10)
  
  ### file 2
  subset2 <- input2[grep(feature, input2$V1),]
  subset2.NoFlanks <- tail(subset2, -10)
  subset2.NoFlanks <- head(subset2.NoFlanks, -10)
  subset2.NoFlanks$coordinates <- 1:nrow(subset2.NoFlanks)   # Make a new column with coordinates starting from 1
  subset2.NoFlanks$Conditions <- condition2  
  subset2.mean <- mean(subset2.NoFlanks$V10)
  
  ### get sno/miRNA gene names rather than IDs
  if (length(args)==4) {
    featureRows <- GTF[grep(feature, GTF$V9),]
    featureRows <- featureRows[1,]
    geneName <- as.character(sub(".*gene_name *(.*?) *; gene_source.*", "\\1", featureRows$V9))
  } else {
    geneName <- feature
  }
  
  ### combine dataframes
  new.df <- rbind(subset1.NoFlanks, subset2.NoFlanks)
  
  ### Plot only if the mean is above the min (arg4)
  if(subset1.mean > mean.cutoff | subset2.mean > mean.cutoff){
  p = ggplot(new.df, aes(x=coordinates, y=V10, fill=Conditions), show.legend = FALSE) +
      geom_ribbon(aes(ymin=V10-V11, ymax=V10+V11, colour=Conditions),
                  alpha=0.3) +
      geom_line(size=0.8) +
      xlab('Nucleotide position') +
      ylab('Read coverage') +
      scale_fill_manual(values=c("red", "blue")) +
      scale_color_manual(values=c("red", "blue")) +
      labs(fill="Conditions", title = geneName)
  plot_list[[feature]] = p
  }
}

#pdf("/home/paul/Documents/Pipelines/tirna-pipeline/NewOutput2/Results/Output2.pdf")
pdf(args[3])
for(feature in featuresUnion) {
  print(plot_list[[feature]])
}
dev.off()


