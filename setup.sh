#!/bin/bash

usage() { echo "Usage: $0 -p CPUs
" 1>&2; }

if [ $# -eq 0 ]; then
    echo "No arguments provided"
    CPUs=1
fi

while getopts ":hp:" o; do
    case "${o}" in
		h)
			usage
			exit
			;;
		p)
			CPUs="$OPTARG"
			;;
		*)
            echo "Error in input parameters!"
			usage
			exit 1
            ;;
    esac
done

# Setup for tsRNAsearch

### Download human genome and GTF file
wget -q http://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz ./ &
wget -q http://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz ./ &

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
	sudo apt install cutadapt
fi
if ! [ -x "$(command -v fastqc)" ]; then
	sudo apt install fastqc
fi
if ! [ -x "$(command -v trim_galore)" ]; then
	curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.tar.gz -o trim_galore.tar.gz
	tar xvzf trim_galore.tar.gz
	mv TrimGalore-0.4.5/trim_galore scripts/
	sed -i -e 's/trim_galore\ /scripts\/trim_galore\ /g' tsRNAsearch.sh
	#echo "trim_galore downloaded and uncompressed. Please add trim_galore to your path" >> Setup_instructions.txt
fi

if ! [ -x "$(command -v hisat2)" ]; then
	sudo apt install hisat2 #2.1.0
fi

if ! [ -x "$(command -v featureCounts)" ]; then
	wget -q https://sourceforge.net/projects/subread/files/subread-1.6.3/subread-1.6.3-Linux-x86_64.tar.gz #featureCounts 1.6.3
	tar xvfz subread-1.6.3-Linux-x86_64.tar.gz
	mv subread-1.6.3-Linux-x86_64/bin/featureCounts scripts/
	sed -i -e 's/featureCounts\ /scripts\/featureCounts\ /g' tsRNAsearch.sh
	#echo "Subreads (incl. bin/featureCounts) downloaded and uncompressed. Please add featureCounts to your path" >> Setup_instructions.txt
fi

wait # Wait for genome to download

### Build genome index
gunzip Homo_sapiens.GRCh37.87.gtf.gz  # Gunzip GTF file
mv Homo_sapiens.GRCh37.87.gtf DBs/
hisat2-build -p $CPUs Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz DBs/hisat2_index/Homo_sapiens.GRCh37.dna.primary_assembly

#echo "
#
#PLEASE READ Setup_instructions.txt IF THIS FILE EXISTS
#
#If it does not exist, you should be ready to execute the tiRNA-pipeline
#"


