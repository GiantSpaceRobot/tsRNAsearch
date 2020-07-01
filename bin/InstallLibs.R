#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### Establish the R libraries required for tsRNAsearch execution
### 
###-------------------------------------------------------------------------------

### Get R version
R.version <- version
R.version.string <- R.version$version.string
R.version.string <- as.numeric(substr(strsplit(R.version.string, " ")[[1]][3], 1, 3))

### If old version of R, do old package install steps
if(R.version.string < 3.5){

	source("https://bioconductor.org/biocLite.R")
	BiocInstaller::biocLite(c("ggplot2", "gplots", "ggrepel", "VennDiagram", "DESeq2", "devtools", "EnhancedVolcano", "metap", "plyr", "dplyr", "stringr", "reshape2", "xtable", "rmarkdown", "knitr"))
	#devtools::install_github('kevinblighe/EnhancedVolcano')

} else {

	if(!require(BiocManager)){
	  install.packages("BiocManager")
	}
	BiocManager::install(c("ggplot2", "gplots", "ggrepel", "VennDiagram", "EnhancedVolcano", "metap", "plyr", "dplyr", "stringr", "reshape2", "xtable", "rmarkdown", "knitr", update = FALSE)
	install.packages("devtools")
	devtools::install_version("latticeExtra", version="0.6-28")
	devtools::install_version("Hmisc", version="4.4-1")
	BiocInstall::install("DESeq2", update = FALSE)
#	if(!require(ggplot2)){
#	  install.packages("ggplot2")
#	  library(ggplot2)
#	} 
#
#	if(!require(gplots)){
#	  install.packages("gplots")
#	  library(gplots)
#	} 
#
#	if(!require(ggrepel)){
#	  install.packages("ggrepel")
#	  library(ggrepel)
#	} 
#
#	if(!require(VennDiagram)){
#	  install.packages("VennDiagram")
#	  library(VennDiagram)
#	}
#
#	if(!require(DESeq2)){
#	  install.packages("DESeq2")
#	  library(DESeq2)
#	} 
#	if(!require(DESeq2)){ 
#	  # If DESeq2 still not loaded, try install from github
#	  install.packages("devtools") 
#	  devtools::install_github("mikelove/DESeq2") 
#	}
#
#	if(!require(EnhancedVolcano)){ 
#	  install.packages("devtools") 
#	  devtools::install_github("kevinblighe/EnhancedVolcano") 
#	}
#
#	if(!require(metap)){
#	  install.packages("metap")
#	  library(metap)
#	}
#
#	if(!require(plyr)){
#	  install.packages("plyr")
#	  library(plyr)
#	}
#
#	if(!require(dplyr)){
#	  install.packages("dplyr")
#	  library(dplyr)
#	}
#
#	if(!require(stringr)){
#	  install.packages("stringr")
#	  library(stringr)
#	}
#
#	if(!require(reshape2)){
#	  install.packages("reshape2")
#	  library(reshape2)
#	}
#
#	if(!require(xtable)){
#	  install.packages("xtable")
#	  library(xtable)
#	}
#
#	if(!require(rmarkdown)){
#	  install.packages("rmarkdown")
#	  library(rmarkdown)
#	}
#
#	if(!require(knitr)){
#	  install.packages("knitr")
#	  library(knitr)
#	}
#
}
