#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script calculates statistics for read coverage distributions
### 
###-------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
} 

#### Input file
input <- read.table(args[1])
#input <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/Full-test/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.stddev")

df <- split( input , f = input$V1 )  # Split dataframe based on column 1 elements
#features <- input$V1   # Group features by name
#featuresUnion <- union(features, features) 

my.list = list("feature", 
               "sum.of.difference",
               "mean.of.difference",
               "stddev.of.difference",
               "coef.of.variation.of.difference",
               "sum.of.percent",
               "mean.of.percent",
               "stddev.of.percent",
               "coef.of.variation.of.percent",
               "area.difference",
               "distribution.score",
               #"ks.test_p.value",
               #"chi.square.test_p.value",
               "Fragment retained in:")
write.table(my.list, 
            file = paste0(args[2], ".all-features.txt"),
            quote = FALSE, 
            append = TRUE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(my.list, 
            file = paste0(args[2], ".different-distributions.txt"),
            quote = FALSE, 
            append = TRUE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

for(subset in df) {
#for (feature in featuresUnion){
  #feature <- "ValCAC"
  feature <- as.character(subset[1,1])
  #subset <- input[grep(feature, input$V1),]
  #sum.total <- sum(subset$V2, subset$V3)
  #subset$V2 <- subset$V2/normalise.factor
  #subset$V3 <- subset$V3/normalise.factor
  subset$V5 <- (subset$V2-subset$V3)
  subset$V6 <- do.call(pmax, subset[2:3]) # Get a column with max read coverage per nucleotide
  sum.total <- sum(subset$V6)
  #subset$V6 <- abs(subset$V5)
  if(all(is.na(subset$V4))) {  # If all values for the percent of relative difference ar NA, do not calculate values
    mean1 <- 0
    mean2 <- 0
    mean.coverage <- 0
    sum.of.difference <- "NA"
    mean.of.difference <- "NA"
    stddev.of.difference <- "NA"
    sum.of.percent <- "NA"
    mean.of.percent <- "NA"
    ks.pvalue <- "NA"
    #chi.pvalue <- "NA"
    stddev.of.percent <- "NA"
    coef.of.variation <- "NA"
    coef.of.var.of.difference <- "NA"
    condition <- "NA"
    area.diff <- "NA"
    distribution.score <- "NA"
  } else {
    mean1 <- mean(subset$V2)
    mean2 <- mean(subset$V3)
    mean.coverage <- mean(c(mean1, mean2))  # Get the overall mean coverage of the two conditions
    sum.of.difference <- sum(abs(subset$V5), na.rm = TRUE)
    sum.of.percent <- sum(abs(subset$V4), na.rm = TRUE)
    mean.of.percent <- mean((subset$V4), na.rm = TRUE)
    ks.output <- ks.test(subset$V2, subset$V3)
    ks.pvalue <- ks.output$p.value
    #chi.output <- chisq.test(subset$V2, subset$V3)
    #chi.pvalue <- chi.output$p.value
    stddev.of.percent <- sd(subset$V4, na.rm = TRUE)
    mean.of.difference <- mean(subset$V5)
    stddev.of.difference <- sd(subset$V5)
    coef.of.variation <- ((abs(stddev.of.percent)/mean.of.percent) * 100)
    coef.of.var.of.difference <- ((stddev.of.difference/mean.of.difference) * 100)
    #area.diff <- (stddev.of.difference*stddev.of.percent)/1000
    area.diff <- (sum.of.difference/sum.total)*100
    distribution.score <- area.diff*stddev.of.percent
    #distribution.score <- (area.diff*mean.coverage)/1000
    if (sum.of.percent > 0) {
      condition <- "Condition 1"
    } else {
      condition <- "Condition 2"
    }
  }
  my.list <- list(feature, 
                  sum.of.difference,
                  mean.of.difference,
                  stddev.of.difference,
                  coef.of.var.of.difference,
                  sum.of.percent,
                  mean.of.percent,
                  stddev.of.percent,
                  coef.of.variation,
                  area.diff,
                  distribution.score,
                  #ks.pvalue,
                  #chi.pvalue,
                  condition)
  write.table(my.list, 
              file = paste0(args[2], ".all-features.txt"),
              append = TRUE, 
              quote = FALSE, 
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  if (mean.coverage > 10) {  # Mean coverage must be over 10
    if (distribution.score > 500) { # Distribution score must be over 50 to be plotted
      write.table(my.list, 
                file = paste0(args[2], ".different-distributions.txt"),
                append = TRUE, 
                quote = FALSE, 
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
    }
  }
}



