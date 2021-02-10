#!/usr/bin/env Rscript

## Script name:StackedBarplots.R
##
## Purpose of script:
##
## Author: Dr Paul Donovan
##
## Date Created: 7-10-20
##
## Copyright: MIT license
##
## Email: pauldonovan@rcsi.com

args = commandArgs(trailingOnly=TRUE)

# args 1 = Directory with tRNA-mapping-info*.tsv files
# args 2 = Output prefix

suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

### Read files
path.to.files <- args[1]

##### 
file.names <- dir(path.to.files, pattern ="\\.tsv")
new.file.names <- grep(pattern = "tRNA-mapping", x = file.names, value = TRUE)

cDataAll <- NULL
#for (i in 1:length(file.names)){
for (i in 1:length(new.file.names)){
  full.path <- paste0(path.to.files, new.file.names[i])
  file <- read.table(full.path, header = TRUE)
  new.name <- gsub(".tsv", "", x = gsub(pattern = "_tRNA-mapping-information", "", x = new.file.names[i]))
  file$sample <- new.name
  file <- head(file, -1)  # Remove row with total
  cDataAll <- rbind(cDataAll, file)
}

### Plot for each tRNA type
pdf(file = paste0(args[2], "_samples-per-tRNA-species.pdf"))
ggplot(cDataAll, aes(fill=sample, y=percent.as.tRNA.total, x=tRNA))+
  geom_bar(position="fill", stat="identity")  +
  ylab("Proportion") +
  xlab("tRNA Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Proportion of each sample by tRNA species")
dev.off()

### Plot for each sample
my.colors <- length(unique(cDataAll$tRNA)) # Get total number of colors required
pdf(file = paste0(args[2], "_tRNA-species-per-sample.pdf"))
ggplot(cDataAll, aes(fill=tRNA, y=percent.as.tRNA.total, x=sample))+
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdYlBu"))(my.colors)) +
  ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Proportion of each tRNA species by sample")
  #scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdYlGn"))(my.colors))
dev.off()


