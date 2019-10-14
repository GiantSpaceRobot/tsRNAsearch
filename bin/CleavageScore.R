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
  stop("Error: Not enough command line arguments provided. At least three file names required.")
} 

#### Input file
input1 <- read.table(args[1])
input2 <- read.table(args[2])

if (length(args)==4) {
  GTF <- read.table(args[4], sep = "\t")
} 

df1 <- split( input1 , f = input1$V1 )  # Split dataframe based on column 1 elements
df2 <- split( input2 , f = input2$V1 )  # Split dataframe based on column 1 elements

results.df <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), c("feature",
                                                                 "mean.coverage",
                                                                 "ratio.5prime",
                                                                 "ratio.3prime",
                                                                 "5vs3.ratio.percent", 
                                                                 "cleavage.score.raw",
                                                                 "standard.dev.condition1",
                                                                 "standard.dev.condition2",
                                                                 "penalty (%)",
                                                                 "cleavage.score.penalty",
                                                                 "cleavage.score")) # Initialise empty dataframe with column headers
count <- 1
for(subset1 in df1) {
  subset2 <- df2[[count]]  # Get the corresponding subset from the second file
  count <- count + 1
  subset1.avg <- mean(subset1$V3)
  subset2.avg <- mean(subset2$V3)
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
  if(cond1.fiveprime.avg=="NaN" | cond1.fiveprime.avg==0){
    cond1.fiveprime.avg <- 1
  } 
  if(cond1.threeprime.avg=="NaN" | cond1.threeprime.avg==0){
    cond1.threeprime.avg <- 1
  } 
  if(cond2.fiveprime.avg=="NaN" | cond2.fiveprime.avg==0){
    cond2.fiveprime.avg <- 1
  } 
  if(cond2.threeprime.avg=="NaN" | cond2.threeprime.avg==0){
    cond2.threeprime.avg <- 1
  } 
  #ratio of 5' compared to 3'
  ratio5prime <- cond1.fiveprime.avg/cond2.fiveprime.avg
  ratio3prime <- cond1.threeprime.avg/cond2.threeprime.avg
  if(ratio5prime=="NaN"){
    ratio5vs3 <- 0
  } else if(ratio3prime=="NaN"){
    ratio5vs3 <- 0
  } else if(ratio5prime >= ratio3prime){
    ratio5vs3 <- (ratio5prime/ratio3prime)*100
  } else {
    ratio5vs3 <- (ratio3prime/ratio5prime)*100
  }
  cleavage.score.raw <- mean.coverage*ratio5vs3
  
  ##### Calculating cleavage score penalty (this is penalty is expressed as the relationship between the mean and stdev)
  
  ### Condition 1
  mean1.sum <- sum(mean(subset1$V3))
  std1.sum <- sum(subset1$V4)
  std1.size <- ifelse(is.na(std1.sum/mean1.sum), 80, std1.sum/mean1.sum) 
    # Calculate the relationship between sum of mean and sum of standard deviation
    # If the mean is zero, use standard deviation of 80% as pseudovalue for std1.size
  
  ### Condition 2
  mean2.sum <- sum(mean(subset2$V3))
  std2.sum <- sum(subset2$V4)
  std2.size <- ifelse(is.na(std2.sum/mean2.sum), 80, std2.sum/mean2.sum) 
    # Calculate the relationship between sum of mean and sum of standard deviation
    # If the mean is zero, use standard deviation of 80% as pseudovalue for std1.size
  
  ##### Average both cond1 and cond2 std/mean relationships:
  ### Linear penalty function
  penalty <- mean(c(std1.size, std2.size))/100 # The amount that the distribution score will be penalised 
  ifelse(std1.size > 80 || std2.size > 80, 
         penalty <- 1, 
         penalty <- penalty) # If either condition have a stdev over 80% of mean, increase penalty to max (1)
  relative.penalty <- cleavage.score.raw*penalty
  cleavage.score <- cleavage.score.raw - relative.penalty
    results.df[nrow(results.df) + 1,] = list(feature,
                                           mean.coverage,
                                           ratio5prime,
                                           ratio3prime,
                                           ratio5vs3,
                                           cleavage.score.raw,
                                           std1.size,
                                           std2.size,
                                           penalty*100,  # % penalty applied to raw cleavage score
                                           relative.penalty, # Actual penalty 
                                           cleavage.score)  # Finalised cleavage score
  

}

results.df <- results.df[order(-as.numeric(results.df$`cleavage.score`)),]

# New method for generating cleavage score that doesn't involve multiplying by mean.coverage (this method produces less crazy cleavage scores)
# I still prefer the other method for the results that make it to the top
#results.df <- results.df[order(as.numeric(results.df$`mean.coverage`)),]
#results.df$numbers <- 1:nrow(results.df)
#results.df$metric <- results.df$`5vs3.ratio.percent` * results.df$numbers
#results.df <- results.df[order(-as.numeric(results.df$`metric`)),]

# Remove crud from full dataframe
newdata <- results.df[complete.cases(results.df), ]  # Remove NAs
newdata <- newdata[!grepl("Inf", newdata$`5vs3.ratio.percent`),] # Remove Inf
newdata <- newdata[newdata$`5vs3.ratio.percent` > 130, ] # Get high 5' / 3' ratios
newdata <- newdata[newdata$`penalty (%)` < 80, ] # Remove features with high penalty
newdata <- newdata[ which(newdata$mean.coverage > 1), ] # Mean RPM must be 10 or more across both conditions for each feature 
newdata <- newdata[order(-newdata$cleavage.score),]
newdata$feature <- factor(newdata$feature, levels = newdata$feature[order(newdata$cleavage.score)])

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

# If there are more than 20 features, show top 20
if(nrow(newdata) > 20){
  newdata.subset <- head(newdata, n = 20)
} else {
  newdata.subset <- newdata
}

pdf.width <- nrow(newdata.subset)*0.2 + 3
pdf(file = paste0(args[3], ".high-cleavage-score.pdf"), width = pdf.width, height = 5)
ggplot(data = newdata.subset, mapping = aes(feature, 
                                            newdata.subset$cleavage.score, 
                                            color=newdata.subset$cleavage.score)) +
  geom_point() +
  ggtitle("Cleavage score") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0, size=8)) +
  scale_color_gradient(low="blue", high="red") +
  #scale_y_continuous(trans='log2') +   # Change y axis to log scale
  scale_x_discrete(limits = (levels(newdata.subset$cleavage.score))) +
  coord_flip() +
  labs(colour = "Cleavage\n   score", 
       x = "ncRNA", 
       y = "Cleavage score", 
       subtitle = "Max number of features shown = 20")
dev.off()


