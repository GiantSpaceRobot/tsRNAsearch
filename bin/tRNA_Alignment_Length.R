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

#setwd("/home/paul/Documents/Pipelines/work/4d/287b8181d217ec2232d474ef7f9eac/")

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
} 

get.lengths <- function(split.line){
  ### Get tRNA type
  trna <- split.line[3]
  trna.type <- strsplit(trna, '-')[[1]][2] 
  ### Get longest perfect matching alignment in each SAM line
  CIGAR <- split.line[6]
  CIGAR.split <- unlist(strsplit(CIGAR[[1]], "(?<=[A-Z])", perl = TRUE))
  CIGAR.split.matching <- CIGAR.split[grep("M", CIGAR.split)]
  CIGAR.longest.match <- max(as.integer(gsub("[A-Z]", "", CIGAR.split.matching))) # Get largest integer
  ### Create row for dataframe
  new.row <- c(trna.type, CIGAR.longest.match)
  return(new.row)
}

length.DF <- data.frame() # Create empty DF
con <- read.delim(args[1], sep = "\t", fill = T, row.names = NULL)
#con <- read.delim("SRR1273998_trimmed_accepted_hits_tRNAs_no-header.sam", sep = "\t", fill = T, row.names = NULL)
#my.test <- con[189445:189450,]
length.DF <- data.frame(t(apply(con, 1, get.lengths))) # Apply function to SAM dataframe

### Clean up dataframe 
names(length.DF) <- c("tRNA", "alignment.length")
length.summary.DF <- t(as.data.frame.matrix(table(length.DF))) # Get table of counts for each tRNA
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