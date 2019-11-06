#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### Generate read alignment length plots 
### 
###-------------------------------------------------------------------------------

library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
} 

myPath <- args[1]
file.names <- dir(myPath, pattern =".txt")
cDataAll <- data.frame(my.rownames=(0))
for (i in 1:length(file.names)){
  full.path <- paste0(myPath, file.names[i])
  my.file <- read.table(full.path, header = TRUE)
  my.file$my.rownames <- rownames(my.file)
  cDataAll <- merge(cDataAll, my.file, by="my.rownames", all = TRUE)
}
rownames(cDataAll) <- cDataAll$my.rownames
intermediate.df1 <- cDataAll %>%
  select(-my.rownames) # Remove rownames column
intermediate.df1 <- intermediate.df1[-1,] %>%
  t() %>%
  data.frame()
intermediate.df1[is.na(intermediate.df1)] <- 0
intermediate.df1$tRNAs <- gsub("\\..*", "", rownames(intermediate.df1))
intermediate.df2 <- ddply(intermediate.df1, "tRNAs", numcolwise(sum))
rownames(intermediate.df2) <- intermediate.df2$tRNAs
length.summary.DF <- intermediate.df2 %>%
  select(-tRNAs) %>%
  t() %>%
  data.frame()
rownames(length.summary.DF) <- as.numeric(gsub(".*X", "", rownames(length.summary.DF)))
length.summary.DF <- as.matrix(length.summary.DF)

### Write table
write.table(paste0(args[1], args[2], "_tRNA-read-alignment-lengths.txt"), x = length.summary.DF, quote = F) # Write length summary as output
### Convert counts to proportions
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
pdf(paste0(args[1], args[2], "_tRNA-read-alignment-lengths.pdf"))
ggplot(data = melted.DF, mapping = aes(x = Row.names, y = value)) +
  geom_line() +
  labs(title = paste0("Read Alignment Lengths (", args[2], ")"),
       x = "Length of read alignment",
       y = "Relative count") +
  facet_wrap(facets = vars(variable))
dev.off()

