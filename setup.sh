#!/bin/bash

usage() { echo "Usage: $0 -p CPUs
" 1>&2; }

#if [ $# -eq 0 ]; then
#    echo "No arguments provided. Defaulting to use 1 CPU. 
#	      Please provide parameter '-p #threads' if you wish to use more than 1."
#    CPUs=1
#fi
CPUs=1 # Default CPU number unless overwritten by parameters provided

while getopts ":hg:p:" o; do
    case "${o}" in
		h)
			usage
			exit
			;;
		g)
			g_option="$OPTARG"
			if [ $g_option = "human" ]; then
				genome="human"
			elif [ $g_option = "mouse" ]; then
				genome="mouse"
			elif [ $g_option = "both" ]; then
				genome="both"
			else
				usage
				echo "Error: the acceptable arguments for the -g parameter are 'human', 'mouse', 'both'"
				exit
			fi
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

### Define functions
function human_genome () {
	### Download human genome
	wget -q http://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz ./ &
	wget -q http://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz ./ &
	wait
	echo "Gunzipping human genome files..."
	gunzip Homo_sapiens.GRCh37.87.gtf.gz &  # Gunzip GTF file
	gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz &
	wait
	python bin/GTF_DuplicateRemover.py additional-files/Homo-sapiens_All-ncRNAs.txt Homo_sapiens.GRCh37.87.gtf Homo_sapiens.GRCh37.87.NoDuplicates.gtf
	rm Homo_sapiens.GRCh37.87.gtf &
	mv Homo_sapiens.GRCh37.87.NoDuplicates.gtf DBs/Homo_sapiens.GRCh37.87.gtf
	echo "Building human genome index..."
	hisat2-build -p $CPUs Homo_sapiens.GRCh37.dna.primary_assembly.fa DBs/hisat2_index/Homo_sapiens.GRCh37.dna.primary_assembly
}

function mouse_genome () {
	### Download mouse genome
	wget -q http://ftp.ensembl.org/pub/release-95/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz ./ &
	wget -q http://ftp.ensembl.org/pub/release-95/gtf/mus_musculus/Mus_musculus.GRCm38.95.gtf.gz ./ &
	wait
	echo "Gunzipping mouse genome files..."
	gunzip Mus_musculus.GRCm38.95.gtf.gz &  # Gunzip GTF file
	gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz &
	wait
	mv Mus_musculus.GRCm38.95.gtf DBs/
	echo "Building mouse genome index..."
	hisat2-build -p $CPUs Mus_musculus.GRCm38.dna.primary_assembly.fa DBs/hisat2_index/Mus_musculus.GRCm38.dna.primary_assembly
}

# HISAT2
echo "Looking for hisat2..."
if ! [ -x "$(command -v hisat2)" ]; then
	sudo apt install hisat2 #2.1.0
else
	echo "hisat2 already installed"
fi

### Download genomes
if [ $genome = "human" ]; then
	### Download human genome and GTF file
	echo "Downloading human genome and GTF files..."
	human_genome &
elif [ $genome = "mouse" ]; then
	### Download mouse genome and GTF file
	echo "Downloading mouse genome and GTF files..."
	mouse_genome &
elif [ $genome = "both" ]; then
	### Download human and mouse genomes and GTF files
	echo "Downloading human and mouse genomes and GTF files..."
	human_genome &
	mouse_genome &
elif [ -z $genome ]; then
	### The genome variable is unset
	echo "Please use the -g option with 'human', 'mouse', or 'both' depending on the type of analyses you intend to run"
	exit 1	
fi

# Setup for tsRNAsearch
echo "

tsRNAsearch setup.

This will take approx. 1 hour per genome (downloading and indexing genomes is a slow process)...

"

# curl
echo "Looking for Curl..."
if ! [ -x "$(command -v curl)" ]; then
	sudo apt install curl
else
	echo "Curl already installed"
fi

# python 
echo "Looking for Python..."
if ! [ -x "$(command -v python)" ]; then
	sudo apt install python
else
	echo "Python already installed"
fi

# pip
echo "Looking for Pip..."
if ! [ -x "$(command -v pip)" ]; then
	sudo apt-get install python-pip
else
	echo "Pip already installed"
fi

# python module numpy
pip install numpy

echo "Looking for R and Rscript..."
### Check for libcurl4-openssl-dev
sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev
if ! [ -x "$(command -v Rscript)" ]; then
	sudo apt install r-base  
else
	echo "Rscript already installed"
fi

# Install necessary R libraries
echo "Installing R libraries..."
sudo Rscript bin/InstallLibs.R

echo "Looking for samtools..."
if ! [ -x "$(command -v samtools)" ]; then
	sudo apt install samtools  # Requires htslib 1.7-2
else
	echo "samtools already installed"
fi

# Trim Galore
echo "Looking for Trim_Galore..."
if ! [ -x "$(command -v cutadapt)" ]; then
	sudo apt install cutadapt
fi
if ! [ -x "$(command -v fastqc)" ]; then
	sudo apt install fastqc
fi
#if ! [ -x "$(command -v trim_galore)" ]; then # Check for global trim_galore
#	if [ ! -f bin/trim_galore ]; then # Check for tsRNAsearch trim_galore
		### Download trim_galore and move it to bin Dir. Edit tsRNAsearch to point to new trim_galore location
		#curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.tar.gz -o trim_galore.tar.gz
		#tar xvzf trim_galore.tar.gz
		#mv TrimGalore-0.4.5/trim_galore bin/
#		sed -i -e 's/trim_galore\ /bin\/trim_galore\ /g' tsRNAsearch.sh
#	fi
#else
#	echo "Trim_Galore aleady installed"
#fi

# FeatureCounts
#echo "Looking for featureCounts..."
#if ! [ -x "$(command -v featureCounts)" ]; then
#	if [ ! -f bin/featureCounts ]; then
		### Download featureCounts and move it to bin Dir. Edit tsRNAsearch to point to new featureCounts location
		#wget -q https://sourceforge.net/projects/subread/files/subread-1.6.3/subread-1.6.3-Linux-x86_64.tar.gz ./ #featureCounts 1.6.3
		#tar xvfz subread-1.6.3-Linux-x86_64.tar.gz
		#mv subread-1.6.3-Linux-x86_64/bin/featureCounts bin/
#		sed -i -e 's/featureCounts\ /bin\/featureCounts\ /g' tsRNAsearch.sh
#	fi
#else
#	echo "featureCounts already installed"
#fi

# Create absolute path for bin files
echo "Creating absolute path for tsRNAsearch 'bin' and 'DBs'"
myPath=$(pwd)
sed -i -e "s~ bin~ ${myPath}\/bin~g" tsRNAsearch.sh # using tilde as delimiter here instead of slash as myPath variable contains slashes
sed -i -e "s~ bin~ ${myPath}\/bin~g" tsRNAsearch_DE.sh
sed -i -e "s~DBs~${myPath}\/DBs~g" tsRNAsearch.sh
sed -i -e "s~DBs~${myPath}\/DBs~g" tsRNAsearch_DE.sh
sed -i -e "s~bin\/trim~${myPath}\/bin\/trim~g" tsRNAsearch.sh
sed -i -e "s~bin\/feat~${myPath}\/bin\/feat~g" tsRNAsearch.sh
sed -i -e "s~additional~${myPath}\/additional~g" tsRNAsearch.sh
sed -i -e "s~additional~${myPath}\/additional~g" tsRNAsearch_DE.sh

wait # Wait for things to finish (genome download and database setup)

echo "Done"
