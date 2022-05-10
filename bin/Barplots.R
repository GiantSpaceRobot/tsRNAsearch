#!/usr/bin/env Rscript

## Script name:Barplots.R
##
## Purpose of script: Make barplots
##
## Author: Dr Paul Donovan
##
## Date Created: 29-5-20
##
## Copyright: MIT license
##
## Email: pauldonovan@rcsi.com

args = commandArgs(trailingOnly=TRUE)

# args 1 = Count file
# args 2 = GTF file
# args 3 = Experiment layout file

file.CSV <- read.csv(args[3], header=FALSE)
Condition1 <- as.character(unique(file.CSV$V2)[1])
Condition2 <- as.character(unique(file.CSV$V2)[2])
lvls.df <- as.data.frame(table(file.CSV$V2))
ReplicateNumber1 <- lvls.df[grep(Condition1, lvls.df$Var1),][1,2]
ReplicateNumber2 <- lvls.df[grep(Condition2, lvls.df$Var1),][1,2]

#-----------------#
#    Libraries    #
#-----------------#

# If libraries not installed, install them
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(reshape2))

### Reading in gene mapping file
GTF <- read.table(args[2], sep = "\t") # Read in GTF for name conversion
GTF$NewNames <- paste0(GTF$V1, " (", as.character(sub(".*gene_name *(.*?) *; .*", "\\1", GTF$V9)), ")") # Add new column with gene names
GTF <- GTF[!duplicated(GTF$NewNames),]  # Remove duplicates based on NewNames

### Read in unnormalised count data
cDataAll <- read.table(args[1], 
                        header = TRUE, row.names = 1)

condition1reps <- paste0(Condition1, 1:ReplicateNumber1)
condition2reps <- paste0(Condition2, 1:ReplicateNumber2)
newcolnames <- c(condition1reps, condition2reps) # Create column names using the provided condition names
colnames(cDataAll) <- (newcolnames)
cDataAll[is.na(cDataAll)] <- 0

  ### Get reads per million (RPM)
new.list <- list()
for(i in colnames(cDataAll)){
  scaling.factor <- sum(cDataAll[i])/1e6
  rpm.column <- cDataAll[i]/scaling.factor
  new.list[[i]] <- rpm.column
}
cDataAll <- do.call("cbind", new.list)

  ###
groups <- factor(x=c(rep(Condition1, ReplicateNumber1), rep(Condition2, ReplicateNumber2)), levels=c(Condition1, Condition2))
  
  ### Create barplots for all features identified
  my.levels <- as.character(levels.default(groups)) # Get names of levels
  level1 <- my.levels[1] # get name of level 1
  level2 <- my.levels[2] # get name of level 2
  level1.df <- data.frame(cDataAll) %>% select(c(1:ReplicateNumber1)) # Subset main count matrix to get condition 1
  level2.df <- data.frame(cDataAll) %>% select(c((ReplicateNumber1+1):(ReplicateNumber1+ReplicateNumber2))) # Subset main count matrix to get condition 2
  subtitle <- "Normalised to reads per million (RPM)" # Create subtitle for plot
  barplot.df <- cDataAll

  barplot.df$features <- rownames(barplot.df)
  GTF$ncRNA.type <- gsub(".*\\gene_biotype ", "", GTF$V9) # Make a new column with ncRNA type (remove everything before gene_biotype in col 9)
  GTF$ncRNA.type <- gsub(";.*", "", GTF$ncRNA.type) # Replace the semicolon at the end of string
  ### Match barplot.df ncRNAs with GTF ncRNA types
  for(row in seq_len(nrow(GTF))) {  
    barplot.df$mapped[ substr( barplot.df$features, 0, nchar(as.character(GTF$V1[row]))) == GTF$V1[row] ] <- GTF$ncRNA.type[row]
  }
  barplot.df$mapped <- as.character(barplot.df$mapped) # Convert from factor to character
  ### Assign ncRNA labels for as many ncRNA types as possible
  for(row in seq_len( nrow(barplot.df))) {
    if(is.na(barplot.df$mapped[row])) {   # If no ncRNA type has been assigned yet:
      if(startsWith(barplot.df$features[row], "ENS") == TRUE) {  # If the feature has an 'ENS' ID, it is a gene or unknown feature. Ignore for now.
        barplot.df$mapped[row] <- "other (ncRNA)"
      } else { # If the feature has no ncRNA type yet and does not start with 'ENS', it is a tRNA
        barplot.df$mapped[row] <- "tRNA"
      }
    } 
  } ### End of ncRNA label assignment loop
  
  ### Reformat
  barplot.df$mapped <- as.factor(barplot.df$mapped)
  grouped.df <- barplot.df %>% group_by(mapped)
  grouped.df <- grouped.df %>% select(-c(features))
  
  ### Calculate percent of each feature per sample
  my.new.df <- data.frame(grouped.df)
  summary.df <- aggregate(. ~ mapped, my.new.df, sum) # Sum each column by the ncRNA feature (mapped)
  rownames(summary.df) <- summary.df$mapped
  summary.df <- summary.df %>% select(-c(mapped)) # Convert from RPM to percent
  percent.df <- round(summary.df/10000, 2)
  flipped.percent.df <- data.frame(t(percent.df))
  write.table(x = flipped.percent.df, file = "Percent-Feature-Per-Sample.tsv", quote = FALSE, sep = "\t", col.names=NA)
  
  ### T test for each feature in both conditions
  my.pval.list <- vector("list", ncol(flipped.percent.df))
  for (i in 1:ncol(flipped.percent.df)){
    my.column <- flipped.percent.df[i]
    my.feature <- colnames(my.column)
    my.column.condition1 <- my.column[1:ReplicateNumber1,]  # Get condition1 percentages
    my.column.condition2 <- as.numeric(na.omit(my.column[ReplicateNumber1+1:nrow(my.column),])) # Get condition2 percentages
  if(length(my.column.condition1)<2 || length(my.column.condition2)<2){
    my.t.test.pval <- 1
    my.pval.list[[i]] <- my.t.test.pval
  } else {
    my.t.test <- t.test(my.column.condition1, my.column.condition2)
    my.t.test.pval <- my.t.test$p.value
    my.pval.list[[i]] <- my.t.test.pval
  }
  }
  #names(my.pval.list) <- colnames(flipped.percent.df)
  my.pval.df <- data.frame(matrix(unlist(my.pval.list), nrow=length(my.pval.list), byrow=T))
  rownames(my.pval.df) <- colnames(flipped.percent.df)
  colnames(my.pval.df) <- "p.value"
  my.pval.df$p.adj <- p.adjust(my.pval.df$p.value)
  write.table(x = my.pval.df, file = "Percent-Feature-Per-Sample_T-tests.tsv", quote = FALSE, sep = "\t", col.names=NA)
  
  ### Melt DF for ggplot2
  melted.df <- melt(data = grouped.df, id.vars = "mapped")
  melted.df <- melted.df %>% mutate("Group" = ifelse(grepl(Condition1, variable), Condition1, Condition2))

  ### Prepare for plot
  df.cols <- length(unique(melted.df$mapped)) # Expand color palette
  mycolors <- colorRampPalette(brewer.pal(11, "RdYlBu"))(df.cols) # Expand RdYlBu to # of colours needed based on DF
  new.melted.df <- melted.df %>% group_by(mapped, variable, Group) %>% summarise(value=sum(value)) # Sum values based on groups (make DF simpler)
  myplot <- ggplot(new.melted.df, aes(variable, value, fill = mapped)) +
    geom_col() +
    scale_fill_manual(values = mycolors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(fill = "Feature", x = "", y = "Reads Per Million (RPM)") + 
    ggtitle("Feature Barplot")
  ggsave(width = 7, filename = "Barplot_All-Samples.pdf", plot = myplot) # Save plot using ggplot2 ggsave

