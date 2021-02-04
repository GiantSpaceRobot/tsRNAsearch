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
	  install.packages("BiocManager", repos='http://cran.rstudio.com/')
	}
	#print("Attempting bulk installation of R packages...")
	#BiocManager::install(c("ggplot2", "gplots", "ggrepel", "VennDiagram", "EnhancedVolcano", "metap", "plyr", "dplyr", "stringr", "reshape2", "xtable", "rmarkdown", "knitr"), update = FALSE)
	# For git2r
	if(!require(git2r)){
		print("Installing git2r...")
		BiocManager::install("git2r", update = FALSE)
		library(git2r)
	}
	# For ggplot2
	if(!require(ggplot2)){
		print("Installing ggplot2...")
		BiocManager::install("ggplot2", update = FALSE)
		library(ggplot2)
	}
	# For ggrepel
	if(!require(ggrepel)){
		print("Installing ggrepel...")
		BiocManager::install("ggrepel", update = FALSE)
		library(ggrepel)
	}
	# For VennDiagram
	if(!require(VennDiagram)){
		print("Installing VennDiagram...")
		BiocManager::install("VennDiagram", update = FALSE)
		library(VennDiagram)
	}
	# For EnhancedVolcano
	if(!require(EnhancedVolcano)){
		print("Installing EnhancedVolcano...")
		BiocManager::install("EnhancedVolcano", update = FALSE)
		library(EnhancedVolcano)
	}
	# For tibble
	if(!require(tibble)){
		print("Installing tibble...")
		BiocManager::install("tibble", update = FALSE)
		library(tibble)
	}
	# For plyr
	if(!require(plyr)){
		print("Installing plyr...")
		BiocManager::install("plyr", update = FALSE)
		library(plyr)
	}
	# For dplyr
	if(!require(dplyr)){
		print("Installing dplyr...")
		BiocManager::install("dplyr", update = FALSE)
		library(dplyr)
	}
	# For stringr
	if(!require(stringr)){
		print("Installing stringr...")
		BiocManager::install("stringr", update = FALSE)
		library(stringr)
	}
	# For reshape2
	if(!require(reshape2)){
		print("Installing reshape2...")
		BiocManager::install("reshape2", update = FALSE)
		library(reshape2)
	}
	# For xtable
	if(!require(xtable)){
		print("Installing xtable...")
		BiocManager::install("xtable", update = FALSE)
		library(xtable)
	}
	# For rmarkdown
	if(!require(rmarkdown)){
		print("Installing rmarkdown...")
		BiocManager::install("rmarkdown", update = FALSE)
		library(rmarkdown)
	}
	# For knitr
	if(!require(knitr)){
		print("Installing knitr...")
		BiocManager::install("knitr", update = FALSE)
		library(knitr)
	}
	# For devtools
	if(!require(devtools)){
		print("Installing devtools...")
		BiocManager::install("devtools", update = FALSE)
		library(devtools)
	}
	#print("Installing devtools...")
	#install.packages("devtools")
	# For gplots
	if(!require(caTools)){
		print("Installing gplots and dependency caTools...")
		devtools::install_version("caTools", version="1.17.1.1")
		BiocManager::install("gplots", update = FALSE)
		library(gplots)
	}
	# For metap:
	if(!require(metap)){
		print("Installing metap and dependency mnormt...")
		devtools::install_version("mnormt", version="1.5-5")
		devtools::install_version("metap", version="0.9")
		#BiocManager::install("metap", update = FALSE)
		library(metap)
	}
	#if(!require(annotate)){
	#	devtools::install_version("annotate", version=)
	#	library(annotate)
	#}
	#if(!require(genefilter)){
	#	install.packages("genefilter")
	#	library(genefilter)
	#}
	# For DESeq2
	if(!require(DESeq2)){
		print("Installing DESeq2 and dependencies (latticeExtra and Hmisc)...")
		devtools::install_version("XML", version="3.98-1.12")
		devtools::install_version("latticeExtra", version="0.6-28")
		devtools::install_version("Hmisc", version="4.1-1")
		BiocManager::install("DESeq2", update=FALSE, dependencies=TRUE)
	}
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
