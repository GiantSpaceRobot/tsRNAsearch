#!/usr/bin/env Rscript

###---------------------------------------------------------
### 
### Compare tsRNA counts between two conditions using DESeq2
### 
###---------------------------------------------------------

if(!require(ggplot2)){
  install.packages("ggplot2", dependencies=TRUE, repos='http://cran.rstudio.com/')
  suppressMessages(library(ggplot2))
}
if(!require(dplyr)){
  install.packages("dplyr", dependencies=TRUE, repos='http://cran.rstudio.com/')
  suppressMessages(library(dplyr))
}
if(!require(DESeq2)){
  install.packages("DESeq2", dependencies=TRUE, repos='http://cran.rstudio.com/')
  suppressMessages(library(DESeq2))
}
if(!require(ggrepel)){
  install.packages("ggrepel", dependencies=TRUE, repos='http://cran.rstudio.com/')
  suppressMessages(library(ggrepel))
}
library(gplots)
library(EnhancedVolcano)
#if(!require(EnhancedVolcano)){
#  install.packages("EnhancedVolcano" dependencies=TRUE, repos='http://cran.rstudio.com/')
#  suppressMessages(library(EnhancedVolcano))
#}

args = commandArgs(trailingOnly=TRUE)

layout.file <- read.csv(args[1], header = FALSE)
#layout.file <- read.csv("/home/paul/Documents/Pipelines/tsRNAsearch/additional-files/Example_ExperimentLayout.csv", header = FALSE)
colnames(layout.file) <- c("Sample", "Group")
layout.file$Names= paste0(layout.file$Group, "___", sub('.f.*','', layout.file$Sample)) #substitute the suffix with nothing
layout.file$ShortNames= sub('.f.*','', layout.file$Sample) #substitute the suffix with nothing


input <- as.character(args[2:length(args)])
#input <- as.character(list.files("/home/paul/Documents/Pipelines/NextFlow/tsRNAsearch/MyOutput/", pattern="tsv", full.names=TRUE))
#input <- sample(input) ### Randomise elements in input (for testing)

### Determine the dimensions of the matrix to be created
#mat.colnum = length(input)
input.temp = read.table(input[1])
mat.rownum = nrow(input.temp)
mat.colnum = length(input)
my.matrix <- matrix(, nrow = mat.rownum, ncol = mat.colnum)  # Generate matrix without data

### Combine files into one matrix. Create list of new sample names
my.sample.names <- list()
for(i in 1:length(input)){
  my.file <- input[i]
  newname=basename(my.file)
  newname=gsub(pattern = "_trimmed_accepted_hits_tRNAs.tsv", replacement = "", x = newname)
  my.condition <- layout.file[grep(newname, layout.file$Sample),]  # Get the condition associated with the .tsv data file using the layout file
  my.condition <- as.character(my.condition$Group)[1]
  my.sample.names[[i]] <- paste0(my.condition, "___", newname)
  input.temp = read.table(my.file)
  my.matrix[,i] <- input.temp$V2  # Fill matrix in
}

### Subset matrix. Label rownames and colnames. Convert to dataframe.
rownames(my.matrix) <- input.temp$V1
my.matrix.subset <- my.matrix[rowSums(my.matrix[])>0,]  # Remove all rows that sum to 0
my.dataframe <- my.matrix.subset %>% as.data.frame()
colnames(my.dataframe) <- my.sample.names

### Reorder columns in my.dataframe to match layout file
my.dataframe <- my.dataframe[layout.file$Names]

### Split dataframe by column names
Condition1 <- as.character(head(layout.file$Group, n=1))
Condition2 <- as.character(tail(layout.file$Group, n=1))
condition.1 <- paste0(Condition1, "___")
condition.2 <- paste0(Condition2, "___")
split.dfs <- sapply(c(condition.1, condition.2),
       function(x) my.dataframe[startsWith(names(my.dataframe),x)],
       simplify = FALSE)
condition.1.df <- split.dfs[[1]]
condition.2.df <- split.dfs[[2]]

### Get no. of replicates in each condition
Condition1.replicate.number <- ncol(condition.1.df)
Condition2.replicate.number <- ncol(condition.2.df)

### Principal component analysis
#Tpm PCA plot (not DESeq2)
RPM <- t(t(my.dataframe)/colSums(my.dataframe))*1e6
inlog <- log(RPM)
colLabel <- c(rep("#E41A1C", Condition1.replicate.number), rep("#377EB8", Condition2.replicate.number))
colTy <- c(rep(1, Condition1.replicate.number), rep(2, Condition2.replicate.number))
d <- dist(t(RPM))
fit=cmdscale(d, eig=TRUE, k=2)
x=fit$points[,1]
y=fit$points[,2]
gg.df <- data.frame(cbind("PC1" = x, "PC2" = y))
gg.df$Group <- c(rep(Condition1, Condition1.replicate.number), rep(Condition2, Condition2.replicate.number))
gg.df$ShortNames <- layout.file$ShortNames
pdf(height = 7, width = 10, "PCA.pdf")
ggplot(gg.df, aes(PC1, PC2, color = Group, shape = Group)) + 
  geom_point() + 
  geom_text_repel(label=gg.df$ShortNames, show.legend = F) +
  ggtitle("Principal Component Analysis")
garbage <- dev.off()

### Create a density plot
pdf("Density-Plot.pdf")
max.density <- 0
for(i in 1:ncol(RPM)){
  my.density <- density(inlog[,i])
  if(max(my.density$y) > max.density){
    max.density <- max(my.density$y)
  }
}
inlog.density.ymax <- max.density*1.2 # Get y max and add 20% 
plot(density(inlog[,1]),
     #ylim=c(0,0.2), 
     ylim=c(0,inlog.density.ymax),
     main="Density plot of counts per gene",
     lty=colTy[1],
     xlab="Log of RPM per gene",
     ylab="Density",
     col=colLabel[1])
for(i in 2:ncol(RPM)){
  lines(density(inlog[,i]), lty=colTy[i], col=colLabel[i])
}
legend("topright", legend=layout.file$ShortNames, lty=colTy, col=colLabel)
garbage <- dev.off()

### Differential Expression Analysis
my.dataframe <- my.dataframe + 1 # Add pseudocount to all cells
dds <- DESeqDataSetFromMatrix(countData = my.dataframe, colData = gg.df, formula(~ Group)) # Compare by Group/Condition
dds <- DESeq(dds)
res <- results(dds, cooksCutoff=FALSE)
res <- res[order(res$padj, decreasing=FALSE),]
results.DF <- data.frame(res)

pdf("MA-plot.pdf")
plotMA(res, ylim=c(-2,2))
dev.off()

### Boxplots of top 9 features
pdf(width = 7, height = 10, 
    "Top-Fragments.pdf")
par(mfrow = c(3,3))# 3 x 3 rows and columns of plots
for (i in 1:9){
  my.fragment <- rownames(results.DF)[i]
  plotCounts(dds, gene=my.fragment, intgroup = "Group")
}
dev.off()

rld <- varianceStabilizingTransformation(dds)
#rld <- rlogTransformation(dds)  # Alternative

### Convert the rld transformation into a matrix
rld.matrx <- assay(rld)   
rld.df <- data.frame(rld.matrx)

### Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
sampleDists <- data.frame(sampleDists)
rownames(sampleDists) <- gsub(pattern = ".*___", replacement = "", rownames(sampleDists))
colnames(sampleDists) <- gsub(pattern = ".*___", replacement = "", colnames(sampleDists))
pdf("Distance-Matrix.pdf")
#par(mar=c(6,4,4,5)+0.1) 
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"), 
          #labRow = layout.file$ShortNames, 
          #labCol = layout.file$ShortNames,
          #ColSideColors=mycols[file.names], 
          #RowSideColors=mycols[file.names],
          #cexRow = 0.8,
          #cexCol = 0.8,
          margins=c(12,10),
          srtCol=45,
          main="Sample Distance Matrix")
garbage <- dev.off()

### Rearrange columns and subset
results.DF$features <- rownames(results.DF) # Add column with gene IDs
results.DF <- results.DF %>% arrange(pvalue) %>% relocate(features)
results.subset.DF <- results.DF %>% filter(pvalue < 0.05) 

#results.DF <- subset(results.DF, select = -c(final)) # Remove unneeded columns
#results.DF$gene.name <- ifelse(is.na(results.DF$gene.name), results.DF$features, results.DF$gene.name)

write.csv(results.DF, file=paste0("DESeq2-output_", Condition1, "-vs-", Condition2, ".csv"), quote = F, row.names = F)

up <- (results.DF[!is.na(results.DF$pvalue) & results.DF$pvalue <= 0.05 &
                   results.DF$log2FoldChange >= 0, ])
down <- (results.DF[!is.na(results.DF$pvalue) & results.DF$pvalue <= 0.05 &    
                     results.DF$log2FoldChange <= 0, ]) 

write.table(x = up,
            sep = ",",
            file=paste0("DESeq2-output-upregulated_", Condition1, "-vs-", Condition2, ".csv"), 
            #col.names=NA,
            row.names=F, 
            quote = FALSE)
write.table(x = down,
            sep = ",",
            file=paste0("DESeq2-output-downregulated_", Condition1, "-vs-", Condition2, ".csv"), 
            #col.names=NA,
            row.names=F, 
            quote = FALSE)

### Create Volcano plot using EnhancedVolcano fro Kevin Blighe
volcano <- EnhancedVolcano(results.DF,
                           lab = results.DF$features,
                           x = 'log2FoldChange',
                           y = 'pvalue',
                           xlim = c(min(results.DF$log2FoldChange, na.rm=TRUE),
                                    max(results.DF$log2FoldChange, na.rm=TRUE)),
                           ylim = c(0, max(-log10(results.DF$pvalue), na.rm=TRUE) + 1),
                           pCutoff = 0.05, 
                           FCcutoff = 0.5)#,
                           #drawConnectors = TRUE,
                           #widthConnectors = 0.1,
                           #subtitle = "Cutoff = Uncorrected pvalue < 0.05")
#labSize = 3.0)
ggsave(width = 7, 
       height = 7, 
       filename = paste0("VolcanoPlot_", Condition1, "-vs-", Condition2, "_pvalue.pdf"), 
       plot = volcano) # Save plot using ggplot2 ggsave (error occured using normal R PDF save)

volcano <- EnhancedVolcano(results.DF,
                           lab = results.DF$features,
                           x = 'log2FoldChange',
                           y = 'padj',
                           xlim = c(min(results.DF$log2FoldChange, na.rm=TRUE),
                                    max(results.DF$log2FoldChange, na.rm=TRUE)),
                           ylim = c(0, max(-log10(results.DF$padj), na.rm=TRUE) + 1),
                           pCutoff = 0.05, 
                           FCcutoff = 0.5)#,
                           #drawConnectors = TRUE,
                           #widthConnectors = 0.1,
                           #subtitle = "Cutoff = Corrected pvalue (padj) < 0.05")
#labSize = 3.0)
ggsave(width = 7, 
       height = 7, 
       filename = paste0("VolcanoPlot_", Condition1, "-vs-", Condition2, "_padj.pdf"), 
       plot = volcano) # Save plot using ggplot2 ggsave (error occured using normal R PDF save)



