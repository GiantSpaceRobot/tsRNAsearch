#!/usr/bin/env Rscript

## Script name:
##
## Purpose of script: Normalise counts by TPM
##
## Author: Dr Paul Donovan
##
## Date Created:
##
## Copyright: MIT license
##
## Email: pauldonovandonegal@gmail.com

args = commandArgs(trailingOnly=TRUE)
#setwd("/home/paul/Documents/Pipelines/work/6a/8cca5d6237622b884cb5bdb6a5063f/")

library(dplyr)

# args 1 = Counts
# args 2 = Total mapped counts
# args 3 = Simple name
# args 4 = tRNA gtf
# args 5 = ncRNA gtf

# tpm function (from Michael Love)
tpm.norm.count.matrix <- function(count.mat, feature.len) {
  x <- count.mat/feature.len
  return(t(t(x)*1e6/colSums(x)))
}
tpm.norm.single.count.vector <- function(count.vec, feature.len) {
  x <- count.vec/feature.len
  return(t(t(x)*1e6/sum(x)))
}

rawdepth = read.csv(args[1], sep = "\t", header = F) # Read input
#rawdepth = read.csv("TotalRNA3_trimmed_accepted_hits_tRNAs_raw.depth", sep = "\t", header = F) # Read input

#mapped = read.csv(args[2], sep = "\t", header = T) # Read input
#mapped = read.csv("Reads_Mapped.txt", sep = "\t", header = T) # Read input

trna.gtf = read.csv(args[3], sep = "\t", header = F) # Read input
#trna.gtf = read.csv("mouse_tRNAs_relative.gtf", sep = "\t", header = F) # Read input

ncrna.gtf = read.csv(args[4], sep = "\t", header = F) # Read input
#ncrna.gtf = read.csv("mouse_ncRNAs_relative_cdhit.gtf", sep = "\t", header = F) # Read input

my.df <- data.frame() # Empty dataframe

## Get length of feature by searching for it in GTF
all.features <- as.character(unique(rawdepth$V1))
for (i in 1:length(all.features)){
  feature <- all.features[i]
  if (grepl("ENS", feature, fixed = T)){
    rows.with.feature <- ncrna.gtf[grepl(feature, ncrna.gtf$V1, fixed = T),]
    median.len.of.feature <- median(rows.with.feature$V5)
  } else {
    rows.with.feature <- trna.gtf[grepl(feature, trna.gtf$V1, fixed = T),]
    median.len.of.feature <- median(rows.with.feature$V5)
  }
  new.row <- data.frame(feature, median.len.of.feature)
  my.df <- rbind(my.df, new.row)
}

colnames(rawdepth) <- c("feature", "position", "count")  # rename columnes
my.merge <- merge(x = rawdepth, y = my.df, by = "feature") # Merge dataframes
colnames(my.merge) <- c("feature", "position", "count", "length")
my.merge$tpm <- round(tpm.norm.single.count.vector(my.merge$count, my.merge$length), 2)

# Get TPM dataframe and write to file
tpm.df <- my.merge %>% select(feature, position, tpm)
write.table(tpm.df, file = args[2], sep = "\t", row.names = F, col.names = F, quote = F)



