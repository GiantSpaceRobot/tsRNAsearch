#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script calculates p-values for all ncRNA features
### 
###-------------------------------------------------------------------------------

library(ggplot2)
library(metap)

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
} 

#### Input files
input1 <- read.table(args[1])
input2 <- read.table(args[1])
GTF <- read.table(args[3], sep = "\t")

col.num1 <- ncol(input1) # Get no. of columns 
input1 <- input1[,-((col.num1 - 1):col.num1), drop = FALSE] #Drop last 2 columns (mean and std)

col.num2 <- ncol(input2) # Get no. of columns 
input2 <- input2[,-((col.num2 - 1):col.num2), drop = FALSE] #Drop last 2 columns (mean and std)
### Read in GTF for name conversions
#if (length(args)==3) {
#  GTF <- read.table(args[3], sep = "\t")
#} 

df1 <- split( input1 , f = input1$V1 )  # Split dataframe based on column 1 elements
df2 <- split( input2 , f = input2$V1 )  # Split dataframe based on column 1 elements

### Initialise empty dataframe with column headers
results.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("feature",
                                                                 "p.value",
                                                                 "p.adj"))  

#test1 <- input1
#test1 <- cbind(subset1, subset2)

### Loop over feature subsets
for(subset1 in df1) {
  feature <- as.character(subset1[1,1])
  subset2 <- df2[[feature]]
  if(startsWith(feature, "ENS")) {
    featureRows <- GTF[grep(feature, GTF$V9),]
    featureRows <- featureRows[1,]
    geneName <- as.character(sub(".*gene_name *(.*?) *; gene_source.*", "\\1", featureRows$V9))
    feature <- paste0(feature," (",geneName,")")
  } 
  rpm1 <- subset(subset1, select = c(rep(FALSE,2), TRUE))
  rpm2 <- subset(subset2, select = c(rep(FALSE,2), TRUE))
  raw.p.vals <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("feature", # Create empty dataframe
                                                                   "p.value"))  
  for(rowNum in rownames(subset1)) { # For each specific row in the feature subset DFs
    condition1 <- as.numeric(rpm1[rowNum,])
    condition2 <- as.numeric(rpm2[rowNum,])
    ### If the number of replicates is greater than 1 for each condition, run t-test, otherwise assign p-value as 1.
    if (length(condition1) == 1 || length(condition2) == 1) { # 
      # Use an alternative to t-test? Or calculate nothing? Go with nothing for now 
      new.p.val <- 1
      raw.p.vals[nrow(raw.p.vals) + 1,] = list(feature, new.p.val)
    } else {
      ### t.test
      t.test.results <- t.test(condition1, condition2, paired = F) # Run a t-test comparing conditions. These data points are not paired
      new.p.val <- ifelse(is.nan(t.test.results$p.value), 1, t.test.results$p.value) # If the p-value couldn't be calculated, use 1 as p-val
      raw.p.vals[nrow(raw.p.vals) + 1,] = list(feature, new.p.val)
      ### Wilcoxon test
      #wilcoxon.results <- wilcox.test(condition1, condition2, paired = T)
      #raw.p.vals[nrow(raw.p.vals) + 1,] = list(feature, wilcoxon.results$p.value)
      #new.p.val <- ifelse(is.nan(wilcoxon.results$p.value), 1, wilcoxon.results$p.value) # If the p-value couldn't be calculated, use 1 as p-val
      #raw.p.vals[nrow(raw.p.vals) + 1,] = list(feature, new.p.val)
    }
  }
  numeric.raw.p.vals <- raw.p.vals[complete.cases(raw.p.vals), ]
  #all.test <- allmetap(numeric.raw.p.vals$p.value, method = "all")
  fishers.method1 <- sumlog(numeric.raw.p.vals$p.value) # Combine p-values using Fisher's method
  #edgington.method1 <- sump(numeric.raw.p.vals$p.value) # Combine p-values using Edgington's method
  new.row <- cbind(feature, "Fishers.method.pvalue" = fishers.method1$p) # Create new row for results dataframe
  results.df <- rbind(results.df, new.row) # Add new to existing dataframe
  
}

write.table(results.df, 
            file = paste0(args[4], "_FisherMethod_pvalues.tsv"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

results.df.copy <- results.df
results.df <- results.df.copy
results.df$Fishers.method.pvalue <- as.numeric(levels(results.df$Fishers.method.pvalue))[results.df$Fishers.method.pvalue]
results.df$negLog10.pval <- as.numeric(format(-log10(results.df$Fishers.method.pvalue), scientific = F, digits = 2))
results.df <- results.df[order(-results.df$negLog10.pval),]
newdata <- results.df[complete.cases(results.df), ]  # Remove NAs
newdata <- newdata[!grepl("Inf", newdata$negLog10.pval),] # Remove Inf
newdata$feature <- factor(newdata$feature, levels = newdata$feature[order(newdata$negLog10.pval)]) # Refactorise to rank order for plotting 

# If there are more than 50 features, show top 20
if(nrow(newdata) > 20){
  newdata.subset <- head(newdata, n = 20)
} else {
  newdata.subset <- newdata
}

write.table(newdata, 
            file = paste0(args[4], ".low.p.values.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

#my.max <- max(newdata.subset$distribution.score)
pdf.width <- nrow(newdata.subset)*0.2 + 3
pdf(file = paste0(args[4], ".low.p.values.pdf"), width = pdf.width, height = 5)
ggplot(data = newdata.subset, mapping = aes(feature, `negLog10.pval`, color=`negLog10.pval`)) +
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

