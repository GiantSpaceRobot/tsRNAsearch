#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script calculates statistics for read coverage distributions
### 
###-------------------------------------------------------------------------------

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
} 

#### Input file
input <- read.table(args[1])

#### Read in GTF for name conversions
if (length(args)==3) {
  GTF <- read.table(args[3], sep = "\t")
} 

df <- split( input , f = input$V1 )  # Split dataframe based on column 1 elements

my.list = list("feature",
               "mean.coverage",
               "sum.of.difference",
               "mean.of.difference",
               "stddev.of.difference",
               "coef.of.variation.of.difference",
               "sum.of.percent",
               "mean.of.percent",
               "stddev.of.percent.relative.difference",
               "coef.of.variation.of.percent",
               "area.difference",
               "distribution.score.raw",
               "distribution.score.penalty",
               "distribution.score",
               "more.reads.in:")
write.table(my.list, 
            file = paste0(args[2], ".distribution-score.all-features.txt"),
            quote = FALSE, 
            append = TRUE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

results.df <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("feature",
                                                                 "mean.coverage",
                                                                 "stddev.of.percent.relative.difference",
                                                                 "area.difference",
                                                                 "distribution.score",
                                                                 "more.reads.in:"))  # Initialise empty dataframe with column headers

for(subset in df) {
  feature <- as.character(subset[1,1])
  ### If feature is a ncRNA (not a tRNA), get gene name
  if(startsWith(feature, "ENS")) {
    featureRows <- GTF[grep(feature, GTF$V9),]
    featureRows <- featureRows[1,]
    geneName <- as.character(sub(".*gene_name *(.*?) *; .*", "\\1", featureRows$V9)) 
    feature <- paste0(feature," (",geneName,")")
  } 
  subset$V7 <- (subset$V2-subset$V3)
  subset$V8 <- do.call(pmax, subset[2:3]) # Get a column with max read coverage per nucleotide
  sum.total <- sum(subset$V8)
  if(all(is.na(subset$V4))) {  # If all values for the percent of relative difference ar NA, do not calculate values
    mean1 <- 0
    mean2 <- 0
    mean.coverage <- 0
    sum.of.difference <- "NA"
    mean.of.difference <- "NA"
    stddev.of.difference <- "NA"
    sum.of.percent <- "NA"
    mean.of.percent <- "NA"
    stddev.of.percent <- "NA"
    coef.of.variation <- "NA"
    coef.of.var.of.difference <- "NA"
    condition <- "NA"
    area.diff <- "NA"
    distribution.score.raw <- "NA"
    relative.penalty <- "NA"
    distribution.score <- "NA"
  } else {
    mean1 <- mean(subset$V2)
    mean2 <- mean(subset$V3)
    mean.coverage <- mean(c(mean1, mean2))  # Get the overall mean coverage of the two conditions
    sum.of.difference <- sum(abs(subset$V7), na.rm = TRUE)
    sum.of.percent <- sum(abs(subset$V4), na.rm = TRUE)
    mean.of.percent <- mean((subset$V4), na.rm = TRUE)
    stddev.of.percent <- sd(subset$V4, na.rm = TRUE)
    if(is.na(stddev.of.percent)) {
      stddev.of.percent = 0
    }
    mean.of.difference <- mean(subset$V7)
    stddev.of.difference <- sd(subset$V7)
    coef.of.variation <- ((abs(stddev.of.percent)/mean.of.percent) * 100)
    coef.of.var.of.difference <- ((stddev.of.difference/mean.of.difference) * 100)
    area.diff <- (sum.of.difference/sum.total)*100
    distribution.score.raw <- area.diff*stddev.of.percent
    
    ##### Calculating distribution score penalty (this is penalty is expressed as the relationship between the mean and stdev)
    
    ### Condition 1
    mean1.sum <- sum(subset$V2)
    std1.sum <- ifelse(sum(subset$V5) > 0, sum(subset$V5), 1) # If the sum of the standard deviation is zero, use 1 as pseudovalue
    std1.size <- ifelse((std1.sum/mean1.sum) == Inf, 0.8, std1.sum/mean1.sum) 
      # Calculate the relationship between sum of mean and sum of standard deviation
      # If the mean is zero, use 0.8 (standard deviation is 80%) as pseudovalue for std1.size

    ### Condition 2
    mean2.sum <- sum(subset$V3)
    std2.sum <- ifelse(sum(subset$V6) > 0, sum(subset$V6), 1) # If the sum of the standard deviation is zero, use 1 as pseudovalue
    std2.size <- ifelse((std2.sum/mean2.sum) == Inf, 0.8, std2.sum/mean2.sum)
      # Calculate the relationship between sum of mean and sum of standard deviation
      # If the mean is zero, use 0.8 (standard deviation is 80%) as pseudovalue for std1.size
    
    ##### Average both cond1 and cond2 std/mean relationships:
    
    ### Linear penalty function
    penalty <- mean(c(std1.size, std2.size)) # The amount that the distribution score will be penalised 
    ifelse(std1.size > 0.9 || std2.size > 0.9, 
           penalty <- 1, 
           penalty <- penalty) # If either condition have a stdev over 90% of mean, increase penalty to max (1)
    relative.penalty <- distribution.score.raw*penalty
    
    ### Exponential penalty function
    #penalty <- (mean(c(std1.size, std2.size)))*10 # The amount that the distribution score will be penalised 
    #penalty <- (penalty^2)/100 
    #relative.penalty <- distribution.score.raw*penalty
    
    ### Distribution score with penalty applied:
    distribution.score <- distribution.score.raw - relative.penalty
    
    ### 
    if (mean.of.difference > 0) {
      condition <- "Condition 1"
    } else {
      condition <- "Condition 2"
    }
  }
  my.list <- list(feature,
                  mean.coverage,
                  sum.of.difference,
                  mean.of.difference,
                  stddev.of.difference,
                  coef.of.var.of.difference,
                  sum.of.percent,
                  mean.of.percent,
                  stddev.of.percent,
                  coef.of.variation,
                  area.diff,
                  distribution.score.raw,
                  relative.penalty,
                  distribution.score,
                  condition)
  write.table(my.list, 
              file = paste0(args[2], ".distribution-score.all-features.txt"),
              append = TRUE, 
              quote = FALSE, 
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  if (mean.coverage > 10) {  # Mean coverage must be over 10
    if (distribution.score > 1000) { # Distribution score must be over 1000 to be plotted
      results.df[nrow(results.df) + 1,] = list(feature,
                                               mean.coverage,
                                               stddev.of.percent,
                                               area.diff,
                                               distribution.score,
                                               condition)
    }
  }
}

### Clean up the results and subset
results.df <- results.df[order(-results.df$distribution.score),]
newdata <- results.df[complete.cases(results.df), ]  # Remove NAs
newdata <- newdata[!grepl("Inf", newdata$distribution.score),] # Remove Inf
newdata$feature <- factor(newdata$feature, levels = newdata$feature[order(newdata$distribution.score)])

# If there are more than 50 features, show top 20
if(nrow(newdata) > 20){
  newdata.subset <- head(newdata, n = 20)
} else {
  newdata.subset <- newdata
}

### Write the high-scoring ncRNAs to file
write.table(newdata, 
            file = paste0(args[2], ".high-distribution-score.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

if (nrow(newdata.subset) < 5){
  pdf.width <- 7
} else {
  pdf.width <- nrow(newdata.subset)*0.2 + 3
}
pdf(file = paste0(args[2], ".high-distribution-score.pdf"), width = pdf.width, height = 5)
ggplot(data = newdata.subset, mapping = aes(feature, `distribution.score`, color=`distribution.score`)) +
  geom_point() +
  #scale_y_continuous(trans='log2') +
  #scale_y_continuous(breaks=seq(0, my.max, by = round(my.max/5))) +
  ggtitle("Distribution Score") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0, size=8)) +
  scale_color_gradient(low="blue", high="red") +
  coord_flip() +
  labs(colour = "Distribution\n    score", 
       x = "ncRNA", 
       y = "Distribution score", 
       subtitle = "Max number of features shown = 20")
dev.off()

