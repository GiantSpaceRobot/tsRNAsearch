#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script calculates p-values for all ncRNA features
### 
###-------------------------------------------------------------------------------

library(ggplot2)
library(metap)
library(stringr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
} 

#### Input files
input1 <- read.table(args[1])
input2 <- read.table(args[2])
GTF <- read.table(args[3], sep = "\t")

col.num1 <- ncol(input1) # Get no. of columns 
input1 <- input1[,-((col.num1 - 1):col.num1), drop = FALSE] #Drop last 2 columns (mean and std)
col.num2 <- ncol(input2) # Get no. of columns 
input2 <- input2[,-((col.num2 - 1):col.num2), drop = FALSE] #Drop last 2 columns (mean and std)

df1 <- split( input1 , f = input1$V1 )  # Split dataframe based on column 1 elements
df2 <- split( input2 , f = input2$V1 )  # Split dataframe based on column 1 elements

### Initialise empty dataframe with column headers
results.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("feature",
                                                                 "p.value",
                                                                 "p.adj"))  
### Initialise empty dataframe with column headers
top.results.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("feature",
                                                                 "p.value",
                                                                 "p.adj"))  

### Loop over feature subsets
for(subset1 in df1) {
  feature <- as.character(subset1[1,1])
  subset2 <- df2[[feature]]
  if(startsWith(feature, "ENS")) {
    featureRows <- GTF[grep(feature, GTF$V9),]
    featureRows <- featureRows[1,]
    geneName <- as.character(sub(".*gene_name *(.*?) *; .*", "\\1", featureRows$V9))
    feature <- paste0(feature," (",geneName,")")
  } 
  rpm1 <- subset(subset1, select = c(rep(FALSE,2), TRUE))
  rpm2 <- subset(subset2, select = c(rep(FALSE,2), TRUE))
  reps.condition1 <- length(as.numeric(rpm1[1,])) #Get no. of replicates in condition 1
  reps.condition2 <- length(as.numeric(rpm2[1,])) #Get no. of replicates in condition 2
  mean1 <- mean(rowSums(rpm1)/reps.condition1)
  mean2 <- mean(rowSums(rpm2)/reps.condition2)
  ### If the number of replicates is greater than 1 for each condition, run t-test, otherwise assign p-value as 1.
  if (reps.condition1 == 1 | reps.condition2 == 1) { 
      # Use an alternative to t-test? Or calculate nothing? Go with nothing for now 
      raw.p.vals <- data.frame("p.value" = matrix(1, nrow = nrow(rpm1), ncol = 1))  # Create column of pval = 1
      raw.p.vals$feature <- feature
    } else {  ### t.test
      ### Compare condition1 and 2 dataframes using t.test with mapply. Convert to DF. Transpose. Convert to DF.  
      mapply.df <- data.frame(t(data.frame(mapply(t.test, data.frame(t(rpm1)), data.frame(t(rpm2)), paired = F, SIMPLIFY = T))))
      pvals <- mapply.df$p.value
      raw.p.vals <- data.frame("p.value" = matrix(unlist(pvals), nrow=length(pvals), byrow=T))
      raw.p.vals$feature <- feature
    }
  numeric.raw.p.vals <- raw.p.vals[complete.cases(raw.p.vals), ] # Remove all Na/NaN/Inf
  if (nrow(numeric.raw.p.vals) == 0) {  # If entire DF was Na/NaN/Inf, use Fisher's method p.value of 1
    new.row <- cbind(feature, "Fishers.method.pvalue" = 1)
  } else {  # If t.test p.values are available, calculate Fisher's method
    #all.test <- allmetap(numeric.raw.p.vals$p.value, method = "all") # Test all p-value combination methods from metap package
    #edgington.method1 <- sump(numeric.raw.p.vals$p.value) # Combine p-values using Edgington's method
    fishers.method1 <- sumlog(numeric.raw.p.vals$p.value) # Combine p-values using Fisher's method
    new.row <- cbind(feature, "Fishers.method.pvalue" = fishers.method1$p) # Create new row for results dataframe
  }
  results.df <- rbind(results.df, new.row) # Add new to existing dataframe
  if (mean1 >= 1 || mean2 >= 1) { # If the mean RPM of condition 1 or 2 is greater than 1, write to DF
    top.results.df <- rbind(top.results.df, new.row) # Add new to existing dataframe
  }
}

write.table(results.df, 
            file = paste0(args[4], "_FisherMethod_pvalues_all-results.tsv"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(top.results.df, 
            file = paste0(args[4], "_FisherMethod_pvalues.tsv"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

top.results.df$Fishers.method.pvalue <- as.numeric(levels(top.results.df$Fishers.method.pvalue))[top.results.df$Fishers.method.pvalue]
top.results.df$negLog10.pval <- as.numeric(format(-log10(top.results.df$Fishers.method.pvalue), scientific = F, digits = 2))
top.results.df <- top.results.df[order(-top.results.df$negLog10.pval),]
newdata <- top.results.df[complete.cases(top.results.df), ]  # Remove NAs
newdata <- newdata[!grepl("Inf", newdata$negLog10.pval),] # Remove Inf
newdata$feature <- factor(newdata$feature, levels = newdata$feature[order(newdata$negLog10.pval)]) # Refactorise to rank order for plotting 
newdata.ncRNAs <- newdata %>% 
  filter(str_detect(feature, "^ENS")) # Subset ncRNAs
newdata.tsRNAs <- newdata %>% 
  filter(!str_detect(feature, "^ENS")) # Subset tsRNAs

### tsRNAs
# If there are more than 20 features, show top 20
if(nrow(newdata.tsRNAs) > 20){
  newdata.tsRNAs.subset <- head(newdata.tsRNAs, n = 20)
} else {
  newdata.tsRNAs.subset <- newdata.tsRNAs
}

write.table(newdata.tsRNAs, 
            file = paste0(args[4], "_tsRNAs.low.p.values.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

pdf.width <- 7
pdf(file = paste0(args[4], "_tsRNAs.low.p.values.pdf"), width = pdf.width, height = 5)
ggplot(data = newdata.tsRNAs.subset, mapping = aes(feature, `negLog10.pval`, color=`negLog10.pval`)) +
  geom_point() +
  #scale_y_continuous(trans='log2') +
  #scale_y_continuous(breaks=seq(0, my.max, by = round(my.max/5))) +
  ggtitle("Combined p-values using Fisher's method") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0, size=8)) +
  scale_color_gradient(low="blue", high="red") +
  coord_flip() +
  labs(colour = "", 
       x = "ncRNA", 
       y = "Fisher's method (-Log10 p-val)", 
       subtitle = "Max number of features shown = 20")
dev.off()

### ncRNAs
# If there are more than 20 features, show top 20
if(nrow(newdata.ncRNAs) > 20){
  newdata.ncRNAs.subset <- head(newdata.ncRNAs, n = 20)
} else {
  newdata.ncRNAs.subset <- newdata.ncRNAs
}

write.table(newdata.ncRNAs, 
            file = paste0(args[4], "_ncRNAs.low.p.values.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

pdf.width <- 7
pdf(file = paste0(args[4], "_ncRNAs.low.p.values.pdf"), width = pdf.width, height = 5)
ggplot(data = newdata.ncRNAs.subset, mapping = aes(feature, `negLog10.pval`, color=`negLog10.pval`)) +
  geom_point() +
  #scale_y_continuous(trans='log2') +
  #scale_y_continuous(breaks=seq(0, my.max, by = round(my.max/5))) +
  ggtitle("Combined p-values using Fisher's method") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0, size=8)) +
  scale_color_gradient(low="blue", high="red") +
  coord_flip() +
  labs(colour = "", 
       x = "ncRNA", 
       y = "Fisher's method (-Log10 p-val)", 
       subtitle = "Max number of features shown = 20")
dev.off()
