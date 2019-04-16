#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
###
### Calculate distribution score penalty based on mean and stdev
###
###-------------------------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args) == 0) {
  stop(
    "Error: Not enough command line arguments provided. Input file and output file names required."
  )
}

#### Input file
input1 <- read.table(args[1])
input2 <- read.table(args[2])

#input1 <- read.table("/home/paul/Documents/Pipelines/tsRNAsearch/Subset_Marion_Ang-vs-Veh/Results/Data/Intermediate-files/Cond1_ENSG00000163597.stdev")
#input2 <- read.table("/home/paul/Documents/Pipelines/tsRNAsearch/Subset_Marion_Ang-vs-Veh/Results/Data/Intermediate-files/Cond2_ENSG00000163597.stdev")

#if (length(args)==4) {
#  GTF <- read.table(args[4], sep = "\t")
#  #GTF <- read.table("/home/paul/Documents/Pipelines/tsRNAsearch/DBs/hg19-snomiRNA_cdhit.gtf", sep = "\t")
#}

df1 <-
  split(input1 , f = input1$V1)  # Split dataframe based on column 1 elements
df2 <-
  split(input2 , f = input2$V1)  # Split dataframe based on column 1 elements

results.df <-
  setNames(data.frame(matrix(ncol = 2, nrow = 0)),
           c("feature",
             "distribution.score.penalty")) # Initialise empty dataframe with column headers
count <- 1
for (subset1 in df1) {
  subset2 <-
    df2[[count]]  # Get the corresponding subset from the second file
  count <- count + 1
  subset1.avg <- mean(subset1$V2)
  subset2.avg <- mean(subset2$V2)
  #subset1.sum <- sum(subset1$V3)
  #subset2.sum <- sum(subset2$V3)
  mean.coverage <- mean(c(subset1.avg, subset2.avg))
  feature <- as.character(subset1[1, 1])
  #subset1.length <- nrow(subset1)
  #half.length <- as.integer(subset1.length/2)
  ### Get gene name for sno/miRNAs
  #if(startsWith(feature, "ENS")) {
  #  featureRows <- GTF[grep(feature, GTF$V9),]
  #  featureRows <- featureRows[1,]
  #  geneName <- as.character(sub(".*gene_name *(.*?) *; gene_source.*", "\\1", featureRows$V9))
  #  feature <- paste0(feature," (",geneName,")")
  #}
  ### Normalise by reads mapped to each condition:
  #subset1$V3 <- subset1$V3/subset1.sum
  #subset2$V3 <- subset2$V3/subset2.sum
  #if (half.length < 150) {
  #  length.penalty <- half.length/150###
  #} else {
  #  length.penalty <- 1
  #}
  ### Condition 1
  mean1.sum <- sum(subset1$V3)
  std1.sum <- sum(subset1$V4)
  std1.size <-
    (std1.sum / mean1.sum) * 100 # Calculate the relationship between sum of mean and sum of standard deviation
  ### Condition 2
  mean2.sum <- sum(subset2$V3)
  std2.sum <- sum(subset2$V4)
  std2.size <-
    (std2.sum / mean2.sum) * 100 # Calculate the relationship between sum of mean and sum of standard deviation
  ### Average both cond1 and cond2 std/mean relationships:
  penalty <- mean(c(std1.size, std2.size))
  results.df[nrow(results.df) + 1,] = list(feature,
                                           penalty)
}

write.table(
  results.df,
  file = args[3],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)
