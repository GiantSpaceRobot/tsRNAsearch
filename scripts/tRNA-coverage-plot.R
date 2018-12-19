#------------------------------------------------------------------------------
# This scripts for making coverage plots from aligned bam files
# You will need the following files to be in the same folder
# Aligned bam and bai files, transcriptome GTF file and exp plan csv file
# 
# You should define the tRNAs of interest in the begining of the script
#
#-------------------------------------------------------------------------------


################################################################################
#                           loading needed libraries                           #
################################################################################
library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)

################################################################################
#                   loading aligned reads and transcriptome                    #
################################################################################

# this is the path to bam,bai, gtf and exp_plan.csv files
path_to_bams <- file.path("/home/paul/Documents/Pipelines/tirna-pipeline")

# add your trna of interest here
trna_list <- c("chr9.tRNA4-ArgTCT", "chr1.tRNA10-AsnGTT", "chr1.tRNA9-ArgTCT")

# load the related exp plan 
# this is used for reading bam files and pooling exp together
sample_table <- read.csv(file.path(path_to_bams, "exp_layout.csv"), header = FALSE)


# finds number of conditions
sample_exp <- unique(sample_table$V1)

### Remove suffix from filenames
sampleNames <- sample_table$V1
sampleNamesShort <- gsub("\\..*","",sampleNames)

# extract genes from transcritome file 
txdb <- makeTxDbFromGFF(file.path(path_to_bams,"DBs/hg19-wholetRNA-CCA.gtf"), format = "gtf", circ_seqs = character())
genes <- exonsBy(txdb, by = "gene")

# create link to bam and bai files
#bam_filenames <- file.path(path_to_bams, paste0("NewOutput/Results/BAMfiles/", sampleNamesShort,"*.bam"))
#bam_filenames <- file.path(path_to_bams, paste0("NewOutput/Results/BAMfiles/", sampleNamesShort, pattern="^.*.bam$"))
#bai_filenames <- file.path(path_to_bams, paste0("NewOutput/Results/BAMfiles/", sampleNamesShort, pattern="^.*bai$"))
bam_filenames <- list.files(path = "NewOutput/Results/BAMfiles/", pattern = "\\.bam$")
bam_filenamesPath <- paste0("NewOutput/Results/BAMfiles/", bam_filenames)
bai_filenames <- list.files(path = "NewOutput/Results/BAMfiles/", pattern = "\\.bai$")
bai_filenamesPath <- paste0("NewOutput/Results/BAMfiles/", bai_filenames)

#bam_filenames <- file.path(path_to_bams, paste0(sample_table$V2,".bam"))
#bai_filenames <- file.path(path_to_bams, paste0(sample_table$V2,".bam.bai"))
bam_files <- BamFileList(bam_filenamesPath, bai_filenamesPath, yieldSize=2000000)

################################################################################
#               find alignment and coverages for tRNAs of interest             #
################################################################################

# only find alignments for tRNAs of interest
param <- ScanBamParam(which = genes[trna_list,])
trna_alignments <- sapply(bam_files, function(x) readGAlignmentsList(x, param = param))

# calcualte coverage for tRNAs of interest
trna_all_coverage <- sapply(trna_alignments, coverage)

################################################################################
#               create coverage plots for each tRNA of interest                #
################################################################################

# loop through tRNAs and plot them based on exp plan groups
for(trna in trna_list){
  trna_coverage <- as.data.frame(sapply(trna_all_coverage, function(x) x[trna,]))
  trna_plot_data <- data.frame(matrix(nrow = nrow(trna_coverage), ncol = length(sample_exp)))
  colnames(trna_plot_data) <- sample_exp
  for(exp in 1:length(sample_exp)){
    trna_plot_data[,exp] <- rowSums(trna_coverage[,grep(paste0("^", sample_exp[exp], ".*value$"), colnames(trna_coverage))])
  }
  trna_gene <- as.data.frame(genes[trna,])
  # cut the introns
  ##trna_plot_data <- trna_plot_data[unlist(sapply(1:nrow(trna_gene), function(x) seq(trna_gene$start[x], trna_gene$end[x]))),]
  # reversing plot data if strand is negative
  if(trna_gene$strand[1] == "-"){
    trna_plot_data <- trna_plot_data[seq(nrow(trna_plot_data), 1, -1),]
  }
  matplot(trna_plot_data, type = "l", lty = 1, ylab = "Coverage", col = c("blue", "red"))
  legend("bottom",names(trna_plot_data), pch = "--", col = c("blue", "red"))
  title(trna_gene$seqnames[1])
  #print(ggplot(reshape2::melt(trna_plot_data), aes(x = rep(1:nrow(trna_plot_data), length(sample_exp)), y = value, color = variable)) + geom_line() +
  #  ggtitle(trna_gene$seqnames[1]) + xlab("") + ylab("Coverage"))
}


