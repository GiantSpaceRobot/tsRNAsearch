#!/bin/bash

usage() { echo "Usage: $0 -p threads
" 1>&2; }

#if [ $# -eq 0 ]; then
#    echo "No arguments provided. Defaulting to use 1 CPU. 
#	      Please provide parameter '-p #threads' if you wish to use more than 1."
#    threads=1
#fi
threads=1 # Default CPU number unless overwritten by parameters provided

while getopts ":hs:p:" o; do
    case "${o}" in
		h)
			usage
			exit
			;;
		s)
			s_option="$OPTARG"
			if [ $s_option = "human" ]; then
				species="human"
			elif [ $s_option = "mouse" ]; then
				species="mouse"
			elif [ $s_option = "both" ]; then
				species="both"
			else
				usage
				echo "Error: the acceptable arguments for the -s parameter are 'human', 'mouse', 'both'"
				exit
			fi
			;;
		p)
			threads="$OPTARG"
			;;
		*)
            echo "Error in input parameters!"
			usage
			exit 1
            ;;
    esac
done

### Define functions
function human_setup () {
	### Download human genome
	wget -q http://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz ./ &
	#wget -q http://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz ./ &
	wait
	echo "Gunzipping human genome files..."
	gunzip Homo_sapiens.GRCh37.87.gtf.gz &  # Gunzip GTF file
	#gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz &
	wait
	python bin/GTF_DuplicateRemover.py additional-files/Homo-sapiens_All-ncRNAs.txt Homo_sapiens.GRCh37.87.gtf Homo_sapiens.GRCh37.87.NoDuplicates.gtf
	rm Homo_sapiens.GRCh37.87.gtf
	mv Homo_sapiens.GRCh37.87.NoDuplicates.gtf DBs/Homo_sapiens.GRCh37.87.gtf
	#echo "Building human genome index..."
	#mkdir -p DBs/genome_index
	#mkdir -p DBs/genome_index/human
	mkdir -p DBs/genome_index/human-ncRNAs
	#STAR --runThreadN $threads --runMode genomeGenerate --genomeDir DBs/genome_index/human/ --genomeFastaFiles Homo_sapiens.GRCh37.dna.primary_assembly.fa 
	bin/STAR --runThreadN $threads --runMode genomeGenerate --genomeDir DBs/genome_index/human-ncRNAs/ --genomeFastaFiles DBs/hg19-combined_tiRNAs_snomiRNAs.fa --genomeSAindexNbases 8
	#echo "Feel free to delete the Homo_sapiens.GRCh37.dna.primary_assembly.fa file in this directory as it is no longer required"
}

function mouse_setup () {
	### Download mouse genome
	#wget -q http://ftp.ensembl.org/pub/release-95/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz ./ &
	wget -q http://ftp.ensembl.org/pub/release-95/gtf/mus_musculus/Mus_musculus.GRCm38.95.gtf.gz ./ &
	wait
	echo "Gunzipping mouse genome files..."
	gunzip Mus_musculus.GRCm38.95.gtf.gz &  # Gunzip GTF file
	#gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz &
	wait
	python bin/GTF_DuplicateRemover.py additional-files/Mus-musculus_All-ncRNAs.txt Mus_musculus.GRCm38.95.gtf Mus_musculus.GRCm38.95.NoDuplicates.gtf
	rm Mus_musculus.GRCm38.95.gtf 
	mv Mus_musculus.GRCm38.95.NoDuplicates.gtf DBs/Mus_musculus.GRCm38.95.gtf
	echo "Building mouse genome index..."
	mkdir -p DBs/genome_index
	#mkdir -p DBs/genome_index/mouse
	mkdir -p DBs/genome_index/mouse-ncRNAs
	#STAR --runThreadN $threads --runMode genomeGenerate --genomeDir DBs/genome_index/mouse/ --genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa
	bin/STAR --runThreadN $threads --runMode genomeGenerate --genomeDir DBs/genome_index/mouse-ncRNAs/ --genomeFastaFiles DBs/GRCm38-combined_tiRNAs_snomiRNAs.fsa
	#echo "Feel free to delete the Mus_musculus.GRCm38.dna.primary_assembly.fa file in this directory as it is no longer required"
}


# STAR
#echo "Looking for STAR..."
#if ! [ -x "$(command -v STAR)" ]; then
#	sudo apt install rna-star
#else
#	echo "STAR already installed"
#fi

### Download genomes
if [ $species = "human" ]; then
	### Download human GTF file
	echo "Downloading human GTF file..."
	human_setup &
elif [ $species = "mouse" ]; then
	### Download mouse GTF file
	echo "Downloading mouse GTF file..."
	mouse_setup &
elif [ $species = "both" ]; then
	### Download human and mouse GTF files
	echo "Downloading human and mouse GTF files..."
	human_setup &
	mouse_setup &
elif [ -z $species ]; then
	### The genome variable is unset
	echo "Please use the -s option with 'human', 'mouse', or 'both' depending on the type of analyses you intend to run"
	exit 1
fi

# Setup for tsRNAsearch
echo "

Beginning tsRNAsearch setup...

"

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

# cutadapt and fastqc
echo "Looking for cutadapt and fastqc (trim_galore)..."
if ! [ -x "$(command -v cutadapt)" ]; then
	sudo apt install cutadapt
fi
if ! [ -x "$(command -v fastqc)" ]; then
	sudo apt install fastqc
fi

# Create absolute path for bin files
echo "Creating absolute path for tsRNAsearch 'bin' and 'DBs'..."
myPath=$(pwd)
sed -i -e "s~ bin~ ${myPath}\/bin~g" bin/tsRNAsearch_single.sh # using tilde as delimiter here instead of slash as myPath variable contains slashes
sed -i -e "s~ bin~ ${myPath}\/bin~g" tsRNAsearch.sh
sed -i -e "s~DBs~${myPath}\/DBs~g" bin/tsRNAsearch_single.sh
sed -i -e "s~DBs~${myPath}\/DBs~g" tsRNAsearch.sh
sed -i -e "s~bin\/trim~${myPath}\/bin\/trim~g" bin/tsRNAsearch_single.sh
sed -i -e "s~bin\/feat~${myPath}\/bin\/feat~g" bin/tsRNAsearch_single.sh
sed -i -e "s~additional~${myPath}\/additional~g" bin/tsRNAsearch_single.sh
sed -i -e "s~additional~${myPath}\/additional~g" tsRNAsearch.sh

wait # Wait for things to finish running

echo "Done"
