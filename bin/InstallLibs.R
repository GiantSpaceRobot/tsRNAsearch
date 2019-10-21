#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### Establish the R libraries required for tsRNAsearch execution
### 
###-------------------------------------------------------------------------------

source("https://bioconductor.org/biocLite.R")

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
} 

if(!require(gplots)){
  install.packages("gplots")
  library(gplots)
} 

if(!require(VennDiagram)){
  install.packages("VennDiagram")
  library(VennDiagram)
}

if(!require(DESeq2)){
  install.packages("DESeq2")
  library(DESeq2)
} 

if(!require(DESeq2)){ 
  # If DESeq2 still not loaded, try install from github
  install.packages("devtools") 
  devtools::install_github("mikelove/DESeq2") 
}

if(!require(EnhancedVolcano)){ 
  install.packages("devtools") 
  devtools::install_github("kevinblighe/EnhancedVolcano") 
}

if(!require(metap)){
  install.packages("metap")
  library(metap)
}

if(!require(plyr)){
  install.packages("plyr")
  library(plyr)
}

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}

if(!require(reshape2)){
  install.packages("reshape2")
  library(reshape2)
}

if(!require(xtable)){
  install.packages("xtable")
  library(xtable)
}
