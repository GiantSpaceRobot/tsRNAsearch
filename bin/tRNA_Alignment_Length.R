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

length.vector <- NULL # Create empty DF
con <- file(args[1], "rt") # Allows file to be read one line at a time
#con <- file("/home/paul/Documents/Pipelines/Analyses_tsRNAsearch/31-10-19_subset_KA2w-vs-KA2w-ctrl/RMK2W_RH_101/tRNA-alignment/tsRNAs_aligned.sam", "rt")
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) { # Read one line at a time (SAM files are usually very large)
  my.string <- oneLine
  if (grepl(pattern = "^@", x = my.string) == T){
  } else {
    my.string
    split.line <- strsplit(my.string, "\t")
    ### Get tRNA type
    trna <- split.line[[1]][3]
    trna.type <- strsplit(trna, '-')[[1]][2] 
    ### Get longest perfect matching alignment in each SAM line
    CIGAR <- split.line[[1]][6]
    CIGAR.split <- unlist(strsplit(CIGAR, "(?<=[A-Z])", perl = TRUE))
    CIGAR.split.matching <- CIGAR.split[grep("M", CIGAR.split)]
    CIGAR.longest.match <- max(as.integer(gsub("[A-Z]", "", CIGAR.split.matching))) # Get largest integer
    ### Add newly created row to dataframe
    new.row <- c(trna.type, CIGAR.longest.match)
    length.vector <- rbind(length.vector, new.row)
  }
} 
close(con)
### Clean up dataframe 
length.DF <- data.frame(length.vector, row.names = NULL)
names(length.DF) <- c("tRNA", "alignment.length")
length.summary.DF <- t(as.data.frame.matrix(table(length.DF))) # Get table of counts for each tRNA
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
pdf(args[2])
ggplot(data = melted.DF, mapping = aes(x = Row.names, y = value)) +
  geom_line() +
  labs(title = "Read Alignment Lengths",
       x = "Length of read alignment",
       y = "Relative count") +
  facet_wrap(facets = vars(variable))
dev.off()