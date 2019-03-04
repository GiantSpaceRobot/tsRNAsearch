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
#input <- read.table("/home/paul/Documents/Pipelines/tsRNAsearch/subset13/Results/Ang_1_ATCACG_L008_R1_001/snomiRNA-alignment/accepted_hits_sorted.depth")

df <- split( input , f = input$V1 )  # Split dataframe based on column 1 elements

results.df <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("feature",
                                                                 "mean.coverage",
                                                                 "mean.fiveprime",
                                                                 "mean.3prime",
                                                                 "5vs3.ratio.percent",
                                                                 "mean.x.ratio"))  # Initialise empty dataframe with column headers

for(subset in df) {
  feature <- as.character(subset[1,1])
  subset.length <- nrow(subset)
  mean.coverage <- mean(subset$V3)
  half.length <- as.integer(subset.length/2)
  #if (half.length < 150) {
  #  length.penalty <- half.length/150###
  #} else {
  #  length.penalty <- 1
  #}
  ### Count no. of zeros in subset. If voer half values are zero, do not calculate ratios.
  zero.percent <- (sum(subset$V3==0)/subset.length)*100
  if(zero.percent >= 50) {
    fiveprime.avg <- 0
    threeprime.avg <- 0
    ratio5vs3 <- 0
    mean.x.ratio <- 0
  } else {
    fiveprime <- subset[1:half.length,]
    five.distribution <- fiveprime$V3
    threeprime <- subset[(half.length+1):subset.length,]
    three.distribution <- threeprime$V3
    fiveprime.avg <- mean(threeprime$V3)
    threeprime.avg <- mean(fiveprime$V3)
    if(fiveprime.avg=="NaN"){
      fiveprime.avg <- 0
    } 
    if(threeprime.avg=="NaN"){
      threeprime.avg <- 0
    } 
    if(fiveprime.avg=="NaN"){
      ratio5vs3 <- 0
    } else if(threeprime.avg=="NaN"){
      ratio5vs3 <- 0
    } else if(fiveprime.avg >= threeprime.avg){
      ratio5vs3 <- (fiveprime.avg/threeprime.avg)*100
    } else {
      ratio5vs3 <- (threeprime.avg/fiveprime.avg)*100
    }
    #if(ratio5vs3=="NaN"){
    #  ratio5vs3 <- 0
    #} 
    mean.x.ratio <- mean.coverage*ratio5vs3
  }
  results.df[nrow(results.df) + 1,] = list(feature,
                                           mean.coverage,
                                           fiveprime.avg,
                                           threeprime.avg,
                                           ratio5vs3,
                                           mean.x.ratio)

}

results.df <- results.df[order(-results.df$mean.x.ratio),]
newdata <- results.df[complete.cases(results.df), ]  # Remove NAs
newdata <- newdata[!grepl("Inf", newdata$mean.x.ratio),] # Remove Inf
newdata <- newdata[newdata$`5vs3.ratio.percent` > 135, ] # Get high 5' / 3' ratios


write.table(results.df, 
            file = paste0(args[2], ".all-features.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
write.table(newdata, 
            file = paste0(args[2], ".potentially-cleaved-features.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

### Make plots
#png.width <- nrow(results.df)*15
#png(file = paste0(args[2], ".all-features.png"), width = png.width, height = 400)
#ggplot(data = results.df, mapping = aes(feature, `-Log10 of KS p-value`, color=`-Log10 of KS p-value`)) +
#  geom_point() +
#  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#  scale_color_gradient(low="blue", high="red")
#dev.off()

pdf.width <- nrow(newdata)*0.2 + 3
pdf(file = paste0(args[2], ".potentially-cleaved-features.pdf"), width = pdf.width, height = 5)
ggplot(data = newdata, mapping = aes(feature, `mean.x.ratio`, color=`mean.x.ratio`)) +
  geom_point() +
  scale_y_continuous(trans='log2') +
  ggtitle("Feature cleavage likelihood") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_gradient(low="blue", high="red") +
  labs(colour = "Cleavage\nscore", 
       x = "ncRNA/gene", 
       y = "Cleavage score", 
       subtitle = "Cleavage score = 5' to 3' ratio (%) multiplied by mean coverage (RPM)")
dev.off()


