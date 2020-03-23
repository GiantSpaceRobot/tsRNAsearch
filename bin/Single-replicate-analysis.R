#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script calculates statistics for 5' vs 3' read coverage distributions
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

results.df <- setNames(data.frame(matrix(ncol = 10, nrow = 0)), c("feature",
                                                                 "mean.coverage",
                                                                 "standard.dev",
                                                                 "mean.fiveprime",
                                                                 "mean.3prime",
                                                                 "5vs3.ratio.percent",
                                                                 "cleavage.score",
                                                                 "distribution.score",
                                                                 "slope.score",
                                                                 "note"))  # Initialise empty dataframe with column headers

for(subset in df) {
  feature <- as.character(subset[1,1])
  subset.length <- nrow(subset)
  ### Slope score calculation
  tRNA.slope.score <- 0
  my.slope.list <- c()
  for(my.row in 1:(nrow(subset) -1)){
    nucleotide1 <- subset[my.row, 3] # Get no. of reads mapped to nucleotide
    nucleotide2 <- subset[my.row +1, 3] # Get no. of reads mapped to next nucleotide
    my.slope <- nucleotide1 - nucleotide2
    my.slope.list <- c(my.slope.list, my.slope)
  }
  tRNA.slope.score <- abs(sum(my.slope.list))
  ### End of slope score calculation
  mean.coverage <- mean(subset$V3)
  subset.std <- sd(subset$V3)
  distribution.score <- mean.coverage*subset.std
  half.length <- as.integer(subset.length/2)
  note <- ""
  if (length(args)==3) {
    featureRows <- GTF[grep(feature, GTF$V9),]
    featureRows <- featureRows[1,]
    geneName <- as.character(sub(".*gene_name *(.*?) *; .*", "\\1", featureRows$V9))
    feature <- paste0(feature," (",geneName,")")
  } 
  ### Count no. of zeros in subset. If over 3/4 values are zero, do not calculate ratios.
  zero.percent <- (sum(subset$V3==0)/subset.length)*100
  if(zero.percent >= 75) {
    fiveprime.avg <- 0
    threeprime.avg <- 0
    ratio5vs3 <- 0
    cleavage.score <- 0
  } else {
    fiveprime <- subset[1:half.length,]
    five.distribution <- fiveprime$V3
    threeprime <- subset[(half.length+1):subset.length,]
    three.distribution <- threeprime$V3
    fiveprime.avg <- mean(fiveprime$V3)
    threeprime.avg <- mean(threeprime$V3)
    ### Convert NaNs into zeros
    if(fiveprime.avg=="NaN"){
      fiveprime.avg <- 0
    } 
    if(threeprime.avg=="NaN"){
      threeprime.avg <- 0
    } 
    ### Generate a ratio number of zero if it is not possible to generate a real one
    if(fiveprime.avg=="0"){
      if(threeprime.avg=="0"){
        ratio5vs3 <- 0
      } else {
        ### Fiveprime average is 0, but threeprime average is a real number
        ratio5vs3 <- threeprime.avg
        note <- "Cleavage score: Real 5' to 3' ratio could not be calculated. Using 3' mean coverage as proxy."
      }
    } else if(threeprime.avg=="0"){
      if(fiveprime.avg=="0"){
      ratio5vs3 <- 0
      } else {
        ### Threeprime average is 0, but fiveprime average is a real number
        ratio5vs3 <- fiveprime.avg
        note <- "Cleavage score: Real 5' to 3' ratio could not be calculated. Using 5' mean coverage as proxy."
      }
    } else if(fiveprime.avg >= threeprime.avg){
      ratio5vs3 <- (fiveprime.avg/threeprime.avg)*100
    } else {
      ratio5vs3 <- (threeprime.avg/fiveprime.avg)*100
    }
    cleavage.score <- mean.coverage*ratio5vs3
  }
  results.df[nrow(results.df) + 1,] = list(feature,
                                           mean.coverage,
                                           subset.std,
                                           fiveprime.avg,
                                           threeprime.avg,
                                           ratio5vs3,
                                           cleavage.score,
                                           distribution.score,
                                           tRNA.slope.score,
                                           note)

}

###
write.table(results.df, 
            file = paste0(args[2], ".all-features.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

### Cleavage score:
results.df <- results.df[order(-results.df$cleavage.score),]
newdata <- results.df[complete.cases(results.df), ]  # Remove NAs
newdata <- newdata[!grepl("Inf", newdata$cleavage.score),] # Remove Inf
newdata <- newdata[newdata$`5vs3.ratio.percent` > 200, ] # Get high 5' / 3' ratios
newdata <- newdata[newdata$mean.coverage > 10, ] # Get high 5' / 3' ratios
newdata$feature <- factor(newdata$feature, levels = newdata$feature[order(newdata$cleavage.score)])

# If there are more than 20 features, show top 20
if(nrow(newdata) > 20){
  newdata.subset <- head(newdata, n = 20)
} else {
  newdata.subset <- newdata
}
pdf.width <- 7

write.table(newdata, 
            file = paste0(args[2], ".high-cleavage-score.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

#pdf.width <- nrow(newdata.subset)*0.2 + 3
pdf(file = paste0(args[2], ".high-cleavage-score.pdf"), width = pdf.width, height = 5)
ggplot(data = newdata.subset, mapping = aes(feature, `cleavage.score`, color=`cleavage.score`)) +
  geom_point() +
  scale_y_continuous(trans='log2') +
  ggtitle("Feature Cleavage Score") +
  coord_flip() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_gradient(low="blue", high="red") +
  labs(colour = "Cleavage\nscore", 
       x = "ncRNA/gene", 
       y = "Cleavage score", 
       subtitle = "Cleavage score = 5' to 3' ratio (%) multiplied by mean coverage (RPM)\nMax number of features shown is 20")
dev.off()


### Distribution score
results.df <- results.df[order(-results.df$distribution.score),]
newdata <- results.df[complete.cases(results.df), ]  # Remove NAs
newdata <- newdata[!grepl("Inf", newdata$distribution.score),] # Remove Inf
newdata <- newdata[newdata$mean.coverage > 10, ] # Get high 5' / 3' ratios
newdata$feature <- factor(newdata$feature, levels = newdata$feature[order(newdata$distribution.score)])

# If there are more than 50 features, show top 50
if(nrow(newdata) > 20){
  newdata.subset <- head(newdata, n = 20)
} else {
  newdata.subset <- newdata
}

write.table(newdata, 
            file = paste0(args[2], ".high-distribution-score.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

#if (nrow(newdata.subset) < 5) {
#	pdf.width <- 7
#} else {
#	pdf.width <- nrow(newdata.subset)*0.2 + 3
#}
pdf.width <- 7

pdf(file = paste0(args[2], ".high-distribution-score.pdf"), width = pdf.width, height = 5)
ggplot(data = newdata.subset, mapping = aes(feature, `distribution.score`, color=`distribution.score`)) +
  geom_point() +
  scale_y_continuous(trans='log2') +
  ggtitle("Feature Distribution Score") +
  coord_flip() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_gradient(low="blue", high="red") +
  labs(colour = "Distribution\nscore", 
       x = "ncRNA/gene", 
       y = "Distribution score", 
       subtitle = "Distribution score = Mean coverage x Standard Deviation\nMax number of features shown is 20")
dev.off()


### Slope score
results.df <- results.df[order(-results.df$slope.score),]
newdata <- results.df[complete.cases(results.df), ]  # Remove NAs
newdata <- newdata[!grepl("Inf", newdata$slope.score),] # Remove Inf
newdata <- newdata[newdata$mean.coverage > 10, ] # Get high 5' / 3' ratios
newdata$feature <- factor(newdata$feature, levels = newdata$feature[order(newdata$slope.score)])

# If there are more than 50 features, show top 50
if(nrow(newdata) > 20){
  newdata.subset <- head(newdata, n = 20)
} else {
  newdata.subset <- newdata
}

write.table(newdata, 
            file = paste0(args[2], ".high-slope-score.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)


pdf.width <- 7
pdf(file = paste0(args[2], ".high-slope-score.pdf"), width = pdf.width, height = 5)
ggplot(data = newdata.subset, mapping = aes(feature, `slope.score`, color=`slope.score`)) +
  geom_point() +
  scale_y_continuous(trans='log2') +
  ggtitle("Feature Slope Score") +
  coord_flip() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_gradient(low="blue", high="red") +
  labs(colour = "Slope\nscore", 
       x = "ncRNA/gene", 
       y = "Slope score", 
       subtitle = "Slope score = Sum of slopes between nucleotides\nMax number of features shown is 20")
dev.off()
