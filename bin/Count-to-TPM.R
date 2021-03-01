#!/usr/bin/env Rscript

## Script name:CombinedScore.R
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
#setwd("/home/paul/Documents/Pipelines/work/9b/8d1c599c217bdc8a2be73ed8467b72/")

library(dplyr)

# args 1 = Counts
# args 2 = Simple name
# args 3 = tRNA gtf
# args 4 = ncRNA gtf

# tpm function (from Michael Love)
tpm.norm.count.matrix <- function(count.mat, feature.len) {
  x <- count.mat/feature.len
  return(t(t(x)*1e6/colSums(x)))
}
tpm.norm.single.count.vector <- function(count.vec, feature.len) {
  x <- count.vec/feature.len
  return(t(t(x)*1e6/sum(x)))
}

counts = read.csv(args[1], sep = "\t", header = F) # Read input
#counts = read.csv("TotalRNA1_trimmed_all-counts.count", sep = "\t", header = F) # Read input

trna.gtf = read.csv(args[3], sep = "\t", header = F) # Read input
#trna.gtf = read.csv("mouse_tRNAs_relative.gtf", sep = "\t", header = F) # Read input

ncrna.gtf = read.csv(args[4], sep = "\t", header = F) # Read input
#ncrna.gtf = read.csv("mouse_ncRNAs_relative_cdhit.gtf", sep = "\t", header = F) # Read input

my.df <- data.frame() # Empty dataframe

for (i in 1:nrow(counts)){
  my.row <- counts[i,]
  feature <- as.character(my.row[,1])
  feature.count <- as.numeric(my.row[,2])
  if (grepl("ENS", feature, fixed = T)){
    rows.with.feature <- ncrna.gtf[grepl(feature, ncrna.gtf$V1, fixed = T),]
    median.len.of.feature <- median(rows.with.feature$V5)
  } else {
    rows.with.feature <- trna.gtf[grepl(feature, trna.gtf$V1, fixed = T),]
    median.len.of.feature <- median(rows.with.feature$V5)
  }
  new.row <- data.frame(feature, feature.count, median.len.of.feature)
  my.df <- rbind(my.df, new.row)
}

colnames(my.df) <- c("feature", "count", "length")
my.df$tpm <- round(tpm.norm.single.count.vector(my.df$count, my.df$length), 2)

# Get TPM dataframe and write to file
tpm.df <- my.df %>% select(feature, tpm)
write.table(tpm.df, file = paste0(args[2], "_tpm.count"), sep = "\t", row.names = F, col.names = F, quote = F)

# Get raw + TPM dataframe and write to file
combined.df <- my.df %>% select(feature, count, tpm)
write.table(tpm.df, file = paste0(args[2], "_raw-tpm.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)



