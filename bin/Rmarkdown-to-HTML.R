#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### Generate HTML from Rmarkdown
### 
###-------------------------------------------------------------------------------


library(rmarkdown)
library(knitr)

args = commandArgs(trailingOnly=TRUE)

rmarkdown::render(input = args[1], 
                  output_format = "all", 
                  output_file = args[2])


