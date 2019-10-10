#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script calculates p-values for all ncRNA features
### 
###-------------------------------------------------------------------------------

#library(ggplot2)
#install.packages("metap")
library(metap)

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
} 

#### Input file
input1 <- read.table(args[1])
input2 <- read.table(args[1])

input1 <- read.table("/home/paul/Documents/Pipelines/tsRNAsearch/Runs/Test_again/Data/Intermediate-files/Everything.cond1.depth")
input2 <- read.table("/home/paul/Documents/Pipelines/tsRNAsearch/Runs/Test_again/Data/Intermediate-files/Everything.cond2.depth")

col.num1 <- ncol(input1) # Get no. of columns 
input1 <- input1[,-((col.num1 - 1):col.num1), drop = FALSE] #Drop last 2 columns (mean and std)

col.num2 <- ncol(input2) # Get no. of columns 
input2 <- input2[,-((col.num2 - 1):col.num2), drop = FALSE] #Drop last 2 columns (mean and std)
#### Read in GTF for name conversions
#if (length(args)==3) {
#  GTF <- read.table(args[3], sep = "\t")
#} 

#input1$mean <- rowMeans(subset(input1, select = c(rep(FALSE,2), TRUE)), na.rm = TRUE)
#input1$stdev <- apply(input1[c(rep(FALSE,2), TRUE)], 1, sd) # Calculate stdev using every third column

df1 <- split( input1 , f = input1$V1 )  # Split dataframe based on column 1 elements
df2 <- split( input2 , f = input2$V1 )  # Split dataframe based on column 1 elements

#write.table(my.list, 
#            file = paste0(args[2], ".all-features.txt"),
#            quote = FALSE, 
#            append = TRUE,
#            sep = "\t",
#            row.names = FALSE,
#            col.names = FALSE)

### Initialise empty dataframe with column headers
results.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("feature",
                                                                 "p.value",
                                                                 "p.adj"))  

#fishers.method <- function(p.value.list) {  #https://www.biostars.org/p/211307/
#  pchisq((sum(log(p.value.list))*-2), df=length(p.value.list)*2, lower.tail=F)
#}

### Loop over feature subsets
for(subset1 in df1) {
  feature <- as.character(subset1[1,1])
  #feature <- "ENSMUSG00000065281"
  #subset1 <- df1[[feature]]
  subset2 <- df2[[feature]]
  #subset1$p.value <- subset(subset1, select = c(rep(FALSE,2), TRUE))
  rpm1 <- subset(subset1, select = c(rep(FALSE,2), TRUE))
  rpm2 <- subset(subset2, select = c(rep(FALSE,2), TRUE))
  #p.val.vector <- NULL
  raw.p.vals <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("feature",
                                                                   "p.value"))  
  for(rowNum in rownames(subset1)) { # For each specific row in the feature subset DFs
    condition1 <- as.numeric(rpm1[rowNum,])
    condition2 <- as.numeric(rpm2[rowNum,])
    ### t.test
    t.test.results <- t.test(condition1, condition2, paired = T)
    new.p.val <- ifelse(is.nan(t.test.results$p.value), 1, t.test.results$p.value) # If the p-value couldn't be calculated, use 1 as p-val
    raw.p.vals[nrow(raw.p.vals) + 1,] = list(feature, new.p.val)
    ### Wilcoxon test
    #wilcoxon.results <- wilcox.test(condition1, condition2, paired = T)
    #raw.p.vals[nrow(raw.p.vals) + 1,] = list(feature, wilcoxon.results$p.value)
    #new.p.val <- ifelse(is.nan(wilcoxon.results$p.value), 1, wilcoxon.results$p.value) # If the p-value couldn't be calculated, use 1 as p-val
    #raw.p.vals[nrow(raw.p.vals) + 1,] = list(feature, new.p.val)
  }
  numeric.raw.p.vals <- raw.p.vals[complete.cases(raw.p.vals), ]
  #all.test <- allmetap(numeric.raw.p.vals$p.value, method = "all")
  fishers.method1 <- sumlog(numeric.raw.p.vals$p.value) # Combine p-values using Fisher's method
  new.row <- cbind(feature, "feature.length" = nrow(subset1), "Fishers.method.pvalue" = fishers.method1$p)
  results.df <- rbind(results.df, new.row)
  
}

write.table(results.df, 
            file = "/home/paul/T.test.feature-length.tsv",
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

#results.df <- results.df.copy
#results.df$feature.length.2 <- as.numeric(levels(results.df$feature.length))[results.df$feature.length]
#results.df$Fishers.method.pvalue.2 <- as.numeric(levels(results.df$Fishers.method.pvalue))[results.df$Fishers.method.pvalue]
#results.df <- results.df[order(results.df$feature.length),]
#results.df$negLog10.pval <- format(-log10(results.df$Fishers.method.pvalue.2), scientific = F, digits = 2)
#super.df <- results.df[results.df$negLog10.pval > 0,]
#plot(super.df$feature.length.2, super.df$negLog10.pval, 
#     log = "x", 
#     main = "ncRNA feature length vs -log10 of Fisher method p-value", 
#     xlab = "Feature Length", ylab = "-log10 p-value (Fisher method)")
#cor(super.df$feature.length.2, as.numeric(super.df$negLog10.pval))


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

write.table(newdata, 
            file = paste0(args[2], ".low.p.value.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

#my.max <- max(newdata.subset$distribution.score)
pdf.width <- nrow(newdata.subset)*0.2 + 3
pdf(file = paste0(args[2], ".low.p.value.pdf"), width = pdf.width, height = 5)
ggplot(data = newdata.subset, mapping = aes(feature, `Fishers.method.pvalue`, color=`Fishers.method.pvalue`)) +
  geom_point() +
  #scale_y_continuous(trans='log2') +
  #scale_y_continuous(breaks=seq(0, my.max, by = round(my.max/5))) +
  ggtitle("Combined p-values using Fishers method") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0, size=8)) +
  scale_color_gradient(low="blue", high="red") +
  coord_flip() +
  labs(colour = "Fishers\nmethod\np-value", 
       x = "ncRNA", 
       y = "Fishers method p-value", 
       subtitle = "Max number of features shown = 20")
dev.off()

