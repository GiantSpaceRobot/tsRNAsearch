#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### Generate read alignment length plots 
### 
###-------------------------------------------------------------------------------

library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
} 

myPath <- "/home/paul/Documents/Pipelines/tsRNAsearch/my.test/"
file.names <- dir(myPath, pattern =".txt")
cDataAll <- NULL
for (i in 1:length(file.names)){
  print(i)
  full.path <- paste0(myPath, file.names[i])
  file <- read.table(full.path, header = TRUE)
  cDataAll <- merge(cDataAll, file, by=0, all=T)
}

write.table(paste0(args[2], ".txt"), x = length.summary.DF, quote = F) # Write length summary as output for further analysis
proportion.DF <- prop.table(length.summary.DF, margin=2)*100 # Get table of proportions
### Generate empty table with complete rownames
max.val <- max(as.integer(rownames(proportion.DF)))
new.rownames <- rep(15:max.val, 1)
emptyDF <- data.frame(matrix(0, nrow = length(new.rownames), ncol = 1))
colnames(emptyDF) <- "empty"
rownames(emptyDF) <- new.rownames
### Merge dataframes and clean
filled.DF <- merge(proportion.DF, emptyDF, by = 0, all = T) # Merge real DF and empty DF to fill rownames
filled.DF[is.na(filled.DF)] <- 0 # Replace NAs with 0
filled.DF <- subset(filled.DF, select=-c(empty)) # Remove empty column
### Reformat for ggplot2
melted.DF <- melt(filled.DF)
melted.DF$Row.names <- as.numeric(melted.DF$Row.names) 
### Plot
pdf(paste0(args[2], ".pdf"))
ggplot(data = melted.DF, mapping = aes(x = Row.names, y = value)) +
  geom_line() +
  labs(title = "Read Alignment Lengths",
       x = "Length of read alignment",
       y = "Relative count") +
  facet_wrap(facets = vars(variable))
dev.off()