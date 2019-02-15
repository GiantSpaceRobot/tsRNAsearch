#!/bin/bash

# Setup for tiRNA pipeline

#python 2.7.15rc1 #numpy #sys
if ! [ -x "$(command -v python)" ]; then
	sudo apt install python
fi

if ! [ -x "$(command -v Rscript)" ]; then
	sudo apt install r-base  #Rscript 3.4.4 #ggplot2 #gplots
fi

### Install necessary R libraries
sudo Rscript scripts/InstallLibs.R

if ! [ -x "$(command -v samtools)" ]; then
	sudo apt install samtools  #1.7 # Requires htslib 1.7-2
fi

# Trim Galore
if ! [ -x "$(command -v cutadapt)" ]; then
	#trim_galore 0.4.4_dev
	sudo apt install cutadapt
fi
if ! [ -x "$(command -v fastqc)" ]; then
	sudo apt install fastqc
fi
if ! [ -x "$(command -v trim_galore)" ]; then
	curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.tar.gz -o trim_galore.tar.gz
	tar xvzf trim_galore.tar.gz
	echo "trim_galore downloaded and uncompressed. Please add trim_galore to your path" >> Setup_instructions.txt
fi

if ! [ -x "$(command -v hisat2)" ]; then
	sudo apt install hisat2 #2.1.0
fi

if ! [ -x "$(command -v Rscript)" ]; then
	wget https://sourceforge.net/projects/subread/files/subread-1.6.3/subread-1.6.3-Linux-x86_64.tar.gz #featureCounts 1.6.3
	tar xvfz subread-1.6.3-Linux-x86_64.tar.gz
	echo "Subreads (incl. bin/featureCounts) downloaded and uncompressed. Please add featureCounts to your path" >> Setup_instructions.txt
fi

echo "

PLEASE READ Setup_instructions.txt IF THIS FILE EXISTS

If it does not exist, you should be ready to execute the tiRNA-pipeline
"


