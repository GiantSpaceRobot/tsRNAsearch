#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### Establish the libraries required for tiRNA-pipeline execution
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

if(!require(DESeq2)){
  install.packages("DESeq2")
  library(DESeq2)
} 

if(!require(DESeq2)){ 
  # If DESeq2 still not loaded, try install from github
  install.packages("devtools") 
  devtools::install_github("mikelove/DESeq2") 
}
