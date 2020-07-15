#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script calculates the slope score
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
#input1 <- read.table("/home/paul/Documents/Pipelines/Analyses_tsRNAsearch/PancreaticCancer_Tumour-vs-Normal_14-7-20/Data/Intermediate-files/DataTransformations/sorted_Everything_ncRNAs.cond1.depth.mean")
#input2 <- read.table("/home/paul/Documents/Pipelines/Analyses_tsRNAsearch/PancreaticCancer_Tumour-vs-Normal_14-7-20/Data/Intermediate-files/DataTransformations/sorted_Everything_ncRNAs.cond2.depth.mean")

if (length(args)==4) {
  GTF <- read.table(args[4], sep = "\t")
} 

df1 <- split( input1 , f = input1$V1 )  # Split dataframe based on column 1 elements
df2 <- split( input2 , f = input2$V1 )  # Split dataframe based on column 1 elements

results.df <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("feature",
                                                                 "mean.coverage",  
                                                                 "slope.score.raw",
                                                                 "standard.dev.condition1",
                                                                 "standard.dev.condition2",
                                                                 "penalty (%)",
                                                                 "slope.score.penalty",
                                                                 "slope.score")) # Initialise empty dataframe with column headers
count <- 1
for(subset1 in df1) {
  subset2 <- df2[[count]]  # Get the corresponding subset from the second file
  count <- count + 1
  subset1.avg <- mean(subset1$V3)
  subset2.avg <- mean(subset2$V3)
  mean.coverage <- mean(c(subset1.avg, subset2.avg))
  feature <- as.character(subset1[1,1])
  subset1.length <- nrow(subset1)
  ### Get gene names for ncRNAs
  if(startsWith(feature, "ENS")) {
    featureRows <- GTF[grep(feature, GTF$V9),]
    featureRows <- featureRows[1,]
    geneName <- as.character(sub(".*gene_name *(.*?) *; .*", "\\1", featureRows$V9))
    feature <- paste0(feature," (",geneName,")")
  }
  
  ### Slope score
  condition1.list <- c()
  condition2.list <- c()
  for(nucleotide in 1:(subset1.length-1)){
    # Condition 1
    cond1.nuc1 <- subset1[nucleotide, 3]
    cond1.nuc2 <- subset1[nucleotide +1, 3]
    cond1.slope <- cond1.nuc1 - cond1.nuc2
    condition1.list <- c(condition1.list, cond1.slope)
    # Condition 2
    cond2.nuc1 <- subset2[nucleotide, 3]
    cond2.nuc2 <- subset2[nucleotide +1, 3]
    cond2.slope <- cond2.nuc1 - cond2.nuc2
    condition2.list <- c(condition2.list, cond2.slope)
  }
  my.new.df <- data.frame(cbind(feature, as.character(1:(subset1.length-1)), condition1.list, condition2.list))
  condition1.slope <- as.numeric(as.character(my.new.df$condition1.list))
  condition2.slope <- as.numeric(as.character(my.new.df$condition2.list))
  my.new.df$slope.subtraction <- condition1.slope - condition2.slope
  slope.score.raw <- sum(abs(my.new.df$slope.subtraction))
  #print(paste0(feature, ",", slope.score.raw))
  ##### Calculating slope score penalty (this is penalty is expressed as the relationship between the mean and stdev)
  
  ### Condition 1
  mean1.sum <- sum(mean(subset1$V3))
  std1.sum <- sum(subset1$V4)
  std1.size <- ifelse(is.na(std1.sum/mean1.sum), 80, std1.sum/mean1.sum) 
  #std1.size <- ifelse(is.na(std1.sum/mean.coverage), 80, std1.sum/mean.coverage) # Alternative penalty based on mean cov
    # Calculate the relationship between sum of mean and sum of standard deviation
    # If the mean is zero, use standard deviation of 80% as pseudovalue for std1.size
  
  ### Condition 2
  mean2.sum <- sum(mean(subset2$V3))
  std2.sum <- sum(subset2$V4)
  std2.size <- ifelse(is.na(std2.sum/mean2.sum), 80, std2.sum/mean2.sum) 
  #std2.size <- ifelse(is.na(std2.sum/mean.coverage), 80, std2.sum/mean.coverage) # Alternative penalty based on mean cov
    # Calculate the relationship between sum of mean and sum of standard deviation
    # If the mean is zero, use standard deviation of 80% as pseudovalue for std1.size
  
  ##### Average both cond1 and cond2 std/mean relationships:
  ### Linear penalty function
  penalty <- mean(c(std1.size, std2.size))/100 # The amount that the distribution score will be penalised 
  #ifelse(std1.size > 80 || std2.size > 80, 
  #       penalty <- 1, 
  #       penalty <- penalty) # If either condition have a stdev over 80% of mean, increase penalty to max (1)
  relative.penalty <- slope.score.raw*penalty
  slope.score <- slope.score.raw - relative.penalty
    results.df[nrow(results.df) + 1,] = list(feature,
                                           mean.coverage,
                                           slope.score.raw,
                                           std1.size,
                                           std2.size,
                                           penalty*100,  # % penalty applied to raw slope score
                                           relative.penalty, # Actual penalty 
                                           slope.score)  # Finalised slope score
  
  
}

results.df <- results.df[order(-as.numeric(results.df$`slope.score`)),]

# Remove crud from full dataframe
newdata <- results.df[complete.cases(results.df), ]  # Remove NAs
newdata <- newdata[!grepl("Inf", newdata$`slope.score`),] # Remove Inf
newdata <- newdata[newdata$`slope.score` > 100, ] # Get high 5' / 3' ratios
newdata <- newdata[newdata$`penalty (%)` < 80, ] # Remove features with high penalty
newdata <- newdata[ which(newdata$mean.coverage > 1), ] # Mean RPM must be 10 or more across both conditions for each feature 
newdata <- newdata[order(-newdata$slope.score),]
newdata$feature <- factor(newdata$feature, levels = newdata$feature[order(newdata$slope.score)])

write.table(results.df, 
            file = paste0(args[3], ".slope-score.all-features.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
write.table(newdata, 
            file = paste0(args[3], ".high-slope-score.txt"),
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

pdf.width <- 7

pdf(file = paste0(args[3], ".high-slope-score.pdf"), width = pdf.width, height = 5)
ggplot(data = newdata.subset, mapping = aes(feature, 
                                            newdata.subset$slope.score, 
                                            color=newdata.subset$slope.score)) +
  geom_point() +
  ggtitle("Slope score") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0, size=8)) +
  scale_color_gradient(low="blue", high="red") +
  #scale_y_continuous(trans='log2') +   # Change y axis to log scale
  scale_x_discrete(limits = (levels(newdata.subset$slope.score))) +
  coord_flip() +
  labs(colour = "Slope\n   score", 
       x = "ncRNA", 
       y = "Slope score", 
       subtitle = "Max number of features shown = 20")
dev.off()


