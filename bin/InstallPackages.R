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

	# For BiocManager
	if(!require(BiocManager)){
	
	  #install.packages("BiocManager")
	  install.packages("https://cran.r-project.org/src/contrib/Archive/BiocManager/BiocManager_1.30.1.tar.gz", repos=NULL)
	}
	# For devtools
	if(!require(devtools)){
		print("Installing devtools...")
		BiocManager::install("devtools", update = FALSE)
		#install.packages("https://cran.r-project.org/src/contrib/Archive/devtools/devtools_2.3.0.tar.gz", repos=NULL)
		library(devtools)
	}
	# For git2r
	if(!require(git2r)){
		print("Installing git2r...")
		install.packages("https://cran.r-project.org/src/contrib/Archive/git2r/git2r_0.27.1.tar.gz", repos=NULL)
		#BiocManager::install("git2r", update = FALSE)
		library(git2r)
	}
	# For ggplot2
	if(!require(ggplot2)){
		print("Installing ggplot2...")
		BiocManager::install("ggplot2", update = FALSE)
		#install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.3.2.tar.gz", repos=NULL)
		#install_version("ggplot2", version = "3.3.2", repos="https://cloud.r-project.org")
		library(ggplot2)
	}
	# For ggrepel
	if(!require(ggrepel)){
		print("Installing ggrepel...")
		BiocManager::install("ggrepel", update = FALSE)
		#install.packages("https://cran.r-project.org/src/contrib/Archive/ggrepel/ggrepel_0.8.2.tar.gz", repos=NULL)
		#install_version("ggrepel", version = "0.8.2", repos="https://cloud.r-project.org")
		library(ggrepel)
	}
	# For VennDiagram
	if(!require(VennDiagram)){
		print("Installing VennDiagram...")
		#BiocManager::install("VennDiagram", update = FALSE)
		#install.packages("https://cran.r-project.org/src/contrib/Archive/VennDiagram/VennDiagram_1.6.19.tar.gz", repos=NULL)
		install_version("VennDiagram", version = "1.6.20", repos="https://cloud.r-project.org")
		library(VennDiagram)
	}
	# For EnhancedVolcano
	if(!require(EnhancedVolcano)){
		print("Installing EnhancedVolcano...")
		BiocManager::install("EnhancedVolcano", update = FALSE)
		#install.packages("", repos=NULL) # No CRAN archive
		library(EnhancedVolcano)
	}
	# For tibble
	if(!require(tibble)){
		print("Installing tibble...")
		#BiocManager::install("tibble", update = FALSE)
		#install.packages("https://cran.r-project.org/src/contrib/Archive/tibble/tibble_3.0.3.tar.gz", repos=NULL)
		install_version("tibble", version = "3.0.3", repos="https://cloud.r-project.org")
		library(tibble)
	}
	# For plyr
	if(!require(plyr)){
		print("Installing plyr...")
		#BiocManager::install("plyr", update = FALSE)
		#install.packages("https://cran.r-project.org/src/contrib/Archive/plyr/plyr_1.8.6.tar.gz", repos=NULL) # Not yet in CRAN
		install_version("plyr", version = "1.8.6", repos="https://cloud.r-project.org")
		library(plyr)
	}
	# For dplyr
	if(!require(dplyr)){
		print("Installing dplyr...")
		#BiocManager::install("dplyr", update = FALSE)
		#install.packages("https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_1.0.0.tar.gz", repos=NULL)
		install_version("dplyr", version = "1.0.0", repos="https://cloud.r-project.org")
		library(dplyr)
	}
	# For stringr
	if(!require(stringr)){
		print("Installing stringr...")
		#BiocManager::install("stringr", update = FALSE)
		#install.packages("https://cran.r-project.org/src/contrib/Archive/stringr/stringr_1.4.0.tar.gz", repos=NULL) # Not yet in CRAN
		install_version("stringr", version = "1.4.0", repos="https://cloud.r-project.org")
		library(stringr)
	}
	# For reshape2
	if(!require(reshape2)){
		print("Installing reshape2...")
		#BiocManager::install("reshape2", update = FALSE)
		#install.packages("https://cran.r-project.org/src/contrib/Archive/reshape2/reshape2_1.4.4.tar.gz", repos=NULL) # Not yet in CRAN
		install_version("reshape2", version = "1.4.4", repos="https://cloud.r-project.org")
		library(reshape2)
	}
	# For xtable
	if(!require(xtable)){
		print("Installing xtable...")
		#BiocManager::install("xtable", update = FALSE)
		#install.packages("https://cran.r-project.org/src/contrib/Archive/xtable/xtable_1.8-4.tar.gz", repos=NULL) # Not yet in CRAN
		install_version("xtable", version = "1.8-4", repos="https://cloud.r-project.org")
		library(xtable)
	}
	# For rmarkdown
	if(!require(rmarkdown)){
		print("Installing rmarkdown...")
		#BiocManager::install("rmarkdown", update = FALSE)
		#install.packages("https://cran.r-project.org/src/contrib/Archive/rmarkdown/rmarkdown_2.3.tar.gz", repos=NULL)
		install_version("rmarkdown", version = "2.3", repos="https://cloud.r-project.org")
		library(rmarkdown)
	}
	# For knitr
	if(!require(knitr)){
		print("Installing knitr...")
		#BiocManager::install("knitr", update = FALSE)
		#install.packages("https://cran.r-project.org/src/contrib/Archive/knitr/knitr_1.29.tar.gz", repos=NULL)
		install_version("knitr", version = "1.29", repos="https://cloud.r-project.org")
		library(knitr)
	}
	# For gplots
	if(!require(caTools)){
		print("Installing gplots and dependency caTools...")
		devtools::install_version("caTools", version="1.17.1.1")
		#install.packages("https://cran.r-project.org/src/contrib/Archive/gplots/gplots_3.0.4.tar.gz", repos=NULL)
		install_version("gplots", version = "3.0.4", repos="https://cloud.r-project.org")
		#BiocManager::install("gplots", update = FALSE)
		library(gplots)
	}
	# For metap:
	if(!require(metap)){
		print("Installing metap and dependency mnormt...")
		devtools::install_version("mnormt", version="1.5-5")
		devtools::install_version("metap", version="0.9")
		BiocManager::install("metap", update = FALSE)
		library(metap)
	}
	# For DESeq2
	if(!require(DESeq2)){
		print("Installing DESeq2 and dependencies (latticeExtra and Hmisc)...")
		devtools::install_version("XML", version="3.98-1.12")
		devtools::install_version("latticeExtra", version="0.6-28")
		devtools::install_version("Hmisc", version="4.1-1")
		install.packages("remotes")
		install_github("mikelove/DESeq2", ref="573fd03", dependencies=TRUE)
		#BiocManager::install("DESeq2", update=FALSE, dependencies=TRUE)
		#install.packages("", repos=NULL)
	}
}
