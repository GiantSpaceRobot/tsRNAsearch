#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script calculates statistics for 5' vs 3' read coverage distributions. Normalise by 
### 
###-------------------------------------------------------------------------------

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
} 

#### Input file
input1 <- read.table(args[1])
input2 <- read.table(args[2])

#input1 <- read.table("/home/paul/Documents/Pipelines/tsRNAsearch/Full_Ang-vs-Veh/Results/Data/Intermediate-files/Distribution-score/sorted_tiRNA.condition1_concatenated.depth.mean")
#input2 <- read.table("/home/paul/Documents/Pipelines/tsRNAsearch/Full_Ang-vs-Veh/Results/Data/Intermediate-files/Distribution-score/sorted_tiRNA.condition2_concatenated.depth.mean")

if (length(args)==4) {
  GTF <- read.table(args[4], sep = "\t")
  #GTF <- read.table("/home/paul/Documents/Pipelines/tsRNAsearch/DBs/hg19-snomiRNA_cdhit.gtf", sep = "\t")
} 

df1 <- split( input1 , f = input1$V1 )  # Split dataframe based on column 1 elements
df2 <- split( input2 , f = input2$V1 )  # Split dataframe based on column 1 elements

results.df <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("feature",
                                                                 "mean.coverage",
                                                                 "ratio.5prime",
                                                                 "ratio.3prime",
                                                                 "5vs3.ratio.percent",
                                                                 "cleavage.score")) # Initialise empty dataframe with column headers
count <- 1
for(subset1 in df1) {
  subset2 <- df2[[count]]  # Get the corresponding subset from the second file
  count <- count + 1
  subset1.avg <- mean(subset1$V3)
  subset2.avg <- mean(subset2$V3)
  #subset1.sum <- sum(subset1$V3)
  #subset2.sum <- sum(subset2$V3)
  mean.coverage <- mean(c(subset1.avg, subset2.avg))
  feature <- as.character(subset1[1,1])
  subset1.length <- nrow(subset1)
  half.length <- as.integer(subset1.length/2)
  ### Get gene name for sno/miRNAs
  if(startsWith(feature, "ENS")) {
    featureRows <- GTF[grep(feature, GTF$V9),]
    featureRows <- featureRows[1,]
    geneName <- as.character(sub(".*gene_name *(.*?) *; gene_source.*", "\\1", featureRows$V9))
    feature <- paste0(feature," (",geneName,")")
  } 
  ### Normalise by reads mapped to each condition:
  #subset1$V3 <- subset1$V3/subset1.sum
  #subset2$V3 <- subset2$V3/subset2.sum
  #if (half.length < 150) {
  #  length.penalty <- half.length/150###
  #} else {
  #  length.penalty <- 1
  #}
  ### Condition1
  cond1.fiveprime <- subset1[1:half.length,]
  cond1.five.distribution <- cond1.fiveprime$V3
  cond1.threeprime <- subset1[(half.length+1):subset1.length,]
  cond1.three.distribution <- cond1.threeprime$V3
  ### Condition2
  cond2.fiveprime <- subset2[1:half.length,]
  cond2.five.distribution <- cond2.fiveprime$V3
  cond2.threeprime <- subset2[(half.length+1):subset1.length,]
  cond2.three.distribution <- cond2.threeprime$V3
  ### Get averages of each
  cond1.fiveprime.avg <- mean(cond1.five.distribution)
  cond1.threeprime.avg <- mean(cond1.three.distribution) 
  cond2.fiveprime.avg <- mean(cond2.five.distribution) 
  cond2.threeprime.avg <- mean(cond2.three.distribution) 
  if(cond1.fiveprime.avg=="NaN"){
    cond1.fiveprime.avg <- 0
  } 
  if(cond1.threeprime.avg=="NaN"){
    cond1.threeprime.avg <- 0
  } 
  if(cond2.fiveprime.avg=="NaN"){
    cond2.fiveprime.avg <- 0
  } 
  if(cond2.threeprime.avg=="NaN"){
    cond2.threeprime.avg <- 0
  } 
  #ratio of 5' compared to 3'
  ratio5prime <- mean(cond1.five.distribution)/mean(cond2.five.distribution)
  ratio3prime <- mean(cond1.three.distribution)/mean(cond2.three.distribution)
  if(ratio5prime=="NaN"){
    ratio5vs3 <- 0
  } else if(ratio3prime=="NaN"){
    ratio5vs3 <- 0
  } else if(ratio5prime >= ratio3prime){
    ratio5vs3 <- (ratio5prime/ratio3prime)*100
  } else {
    ratio5vs3 <- (ratio3prime/ratio5prime)*100
  }
  cleavage.score <- mean.coverage*ratio5vs3
    results.df[nrow(results.df) + 1,] = list(feature,
                                           mean.coverage,
                                           ratio5prime,
                                           ratio3prime,
                                           ratio5vs3,
                                           cleavage.score)
  #}
  

}

results.df <- results.df[order(-as.numeric(results.df$`5vs3.ratio.percent`)),]

# Remove crud from full dataframe
newdata <- results.df[complete.cases(results.df), ]  # Remove NAs
newdata <- newdata[!grepl("Inf", newdata$`5vs3.ratio.percent`),] # Remove Inf
newdata <- newdata[newdata$`5vs3.ratio.percent` > 135, ] # Get high 5' / 3' ratios
#low.ratio.df <- newdata[newdata$`5vs3.ratio.percent` < 66.667, ] # Get low 5' / 3' ratios
#low.ratio.df <- low.ratio.df[low.ratio.df$`5vs3.ratio.percent` > 0, ] # Remove ratios of 0
#newdata <- rbind(high.ratio.df, low.ratio.df) # Combine high and low dataframes
newdata <- newdata[ which(newdata$mean.coverage > 1), ] # Mean RPM must be 10 or more across both conditions for each feature 
newdata <- newdata[order(-newdata$cleavage.score),]


write.table(results.df, 
            file = paste0(args[3], ".all-features.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
write.table(newdata, 
            file = paste0(args[3], ".high-cleavage-score.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

# If there are more than 50 features, show top 50
if(nrow(newdata) > 20){
  newdata.subset <- head(newdata, n = 20)
} else {
  newdata.subset <- newdata
}

pdf.width <- nrow(newdata)*0.2 + 3
pdf(file = paste0(args[3], ".high-cleavage-score.pdf"), width = pdf.width, height = 5)
#ggplot(data = newdata, mapping = aes(feature, newdata$`5vs3.ratio.percent`, color=newdata$`5vs3.ratio.percent`)) +
#  geom_point() +
#  ggtitle("Feature cleavage likelihood") +
#  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#  scale_color_gradient(low="blue", high="red") +
#  labs(colour = "Cleavage\nlikelihood", 
#       x = "ncRNA/gene", 
#       y = "Difference between 5' and 3' ratios (%)", 
#       subtitle = NULL)

#rownames(newdata.subset) <- seq(length=nrow(newdata.subset))


ggplot(data = newdata.subset, mapping = aes(feature, 
                                            newdata.subset$cleavage.score, 
                                            color=newdata.subset$cleavage.score)) +
  geom_point() +
  ggtitle("Feature cleavage analysis") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7)) +
  scale_color_gradient(low="blue", high="red") +
  #scale_y_continuous(trans='log2') +   # Change y axis to log scale
  scale_x_discrete(limits = (levels(newdata.subset$cleavage.score))) +
  labs(colour = "Cleavage\nscore", 
       x = "ncRNA/gene", 
       y = "Cleavage score", 
       subtitle = "Cleavage score = Ratio between 5' ratio (condition 1 vs 2) and\n3' ratio (condition 1 vs 2) multiplied by mean coverage (RPM)\nMax number of features shown is 20")
dev.off()





