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
				echo "Error: the acceptable arguments for the -g parameterhuman are 'human', 'mouse', 'both'"
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
	mv Homo_sapiens.GRCh37.87.gtf DBs/
	echo "Building human genome index..."
	hisat2-build -p $CPUs Homo_sapiens.GRCh37.dna.primary_assembly.fa DBs/hisat2_index/Homo_sapiens.GRCh37.dna.primary_assembly
}

function mouse_genome () {
	### Download mouse genome
	wget -q ftp://ftp.ensembl.org/pub/release-95/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fasta.gz ./ &
	wget -q ftp://ftp.ensembl.org/pub/release-95/gtf/mus_musculus/Mus_musculus.GRCm38.95.gtf.gz ./ &
	wait
	echo "Gunzipping mouse genome files..."
	gunzip Mus_musculus.GRCm38.95.gtf.gz &  # Gunzip GTF file
	gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz &
	wait
	mv Mus_musculus.GRCm38.95.gtf DBs/
	echo "Building mouse genome index..."
	hisat2-build -p $CPUs Mus_musculus.GRCm38.dna.primary_assembly.fa DBs/hisat2_index/Mus_musculus.GRCm38.dna.primary_assembly
}

# Setup for tsRNAsearch
echo "

tsRNAsearch setup.

This will take approx. 1 hour per genome (downloading and indexing genomes is a slow process)...

"

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
fi

# python 2.7.15rc1 
echo "Looking for Python..."
if ! [ -x "$(command -v python)" ]; then
	sudo apt install python
else
	#PyVersion=$(python -V) # Get python version
	echo "Python already installed"
fi

echo "Looking for Rscript..."
if ! [ -x "$(command -v Rscript)" ]; then
	sudo apt install r-base  #Rscript 3.4.4 
else
	echo "Rscript already installed"
fi

### Install necessary R libraries
echo "Installing R libraries..."
sudo Rscript scripts/InstallLibs.R

echo "Looking for samtools..."
if ! [ -x "$(command -v samtools)" ]; then
	sudo apt install samtools  #1.7 # Requires htslib 1.7-2
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
if ! [ -x "$(command -v trim_galore)" ]; then # Check for global trim_galore
	if [ ! -f scripts/trim_galore ]; then # Check for tsRNAsearch trim_galore
		### Download trim_galore and move it to scripts Dir. Edit tsRNAsearch to point to new trim_galore location
		curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.tar.gz -o trim_galore.tar.gz
		tar xvzf trim_galore.tar.gz
		mv TrimGalore-0.4.5/trim_galore scripts/
		sed -i -e 's/trim_galore\ /scripts\/trim_galore\ /g' tsRNAsearch.sh
	fi
else
	echo "Trim_Galore aleady installed"
fi

echo "Looking for hisat2..."
if ! [ -x "$(command -v hisat2)" ]; then
	sudo apt install hisat2 #2.1.0
else
	echo "hisat2 already installed"
fi

echo "Looking for featureCounts..."
if ! [ -x "$(command -v featureCounts)" ]; then
	if [ ! -f scripts/featureCounts ]; then
		### Download featureCounts and move it to scripts Dir. Edit tsRNAsearch to point to new featureCounts location
		wget -q https://sourceforge.net/projects/subread/files/subread-1.6.3/subread-1.6.3-Linux-x86_64.tar.gz ./ #featureCounts 1.6.3
		tar xvfz subread-1.6.3-Linux-x86_64.tar.gz
		mv subread-1.6.3-Linux-x86_64/bin/featureCounts scripts/
		sed -i -e 's/featureCounts\ /scripts\/featureCounts\ /g' tsRNAsearch.sh
	fi
else
	echo "featureCounts already installed"
fi

wait # Wait for things to finish (genome download and database setup)

echo "Done"
