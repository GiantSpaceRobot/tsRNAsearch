library(gplots)

args = commandArgs(trailingOnly=TRUE)

input1 <- read.table(args[1])
input2 <- read.table(args[2])
input3 <- read.table(args[3])

#input1 <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/subset3/Results/Data/DE_Results/DESeq2/DEGs_names-only.txt")
#input2 <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/subset3/Results/Data/Intermediate-files/Distribution-score/High-distribution-scores_feature-names.txt")
#input3 <- read.table("/home/paul/Documents/Pipelines/tirna-pipeline/subset3/Results/Data/Intermediate-files/Potentially-cleaved-features_feature-names.txt")

pdf(file = args[4], 
    width = 8, 
    height = 8)
venn(list(Differentially.Expressed=input1$V1, 
          High.Distribution=input2$V1,
          Potentially.Cleaved=input3$V1))
dev.off()



