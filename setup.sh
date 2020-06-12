#!/bin/bash

usage() { echo "Usage: $0 -s human/mouse/rat/all -p threads
" 1>&2; }

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
			elif [ $s_option = "rat" ]; then
				species="rat"
			elif [ $s_option = "all" ]; then
				species="all"
			else
				usage
				echo "Error: the acceptable arguments for the -s parameter are 'human', 'mouse', 'rat', 'all'"
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

### If no command line arguments provided, quit
if [ -z "$*" ] ; then
	usage
	echo "Error: no command line parameters provided!"
	exit 1
fi

### Define functions
function human_setup () {
	mkdir -p DBs/species_index/human-ncRNAs
	bin/STAR --runThreadN $threads --runMode genomeGenerate --genomeDir DBs/species_index/human-ncRNAs/ --genomeFastaFiles DBs/human_tRNAs-and-ncRNAs_relative_cdhit.fa --genomeSAindexNbases 8
}
function mouse_setup () {
	mkdir -p DBs/species_index/mouse-ncRNAs
	bin/STAR --runThreadN $threads --runMode genomeGenerate --genomeDir DBs/species_index/mouse-ncRNAs/ --genomeFastaFiles DBs/mouse_tRNAs-and-ncRNAs_relative_cdhit.fa
}
function rat_setup () {
	mkdir -p DBs/species_index/rat-ncRNAs
	bin/STAR --runThreadN $threads --runMode genomeGenerate --genomeDir DBs/species_index/rat-ncRNAs/ --genomeFastaFiles DBs/rat_tRNAs-and-ncRNAs_relative_cdhit.fa
}
### Download species data
mkdir -p DBs/species_index
if [ $species = "human" ]; then
	echo "Setting up human database..."
	human_setup &
elif [ $species = "mouse" ]; then
	echo "Setting up mouse database..."
	mouse_setup &
elif [ $species = "rat" ]; then
	echo "Setting up rat database..."
	rat_setup &
elif [ $species = "all" ]; then
	echo "Setting up all available databases..."
	human_setup &
	mouse_setup &
	rat_setup &
elif [ -z $species ]; then
	### The species variable is unset
	echo "Please use the -s option with 'human', 'mouse', 'rat', or 'all' depending on the type of analyses you intend to run"
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
echo "Looking for numpy..."
python -c "import numpy"
if [ $(echo $?) == 1 ]; then # If python numpy call == 1, numpy is not installed
	pip install numpy
else
	echo "Numpy already installed"
fi


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

# pdfunite
echo "Looking for pdfunite..."
if ! [ -x "$(command -v pdfunite)" ]; then
	sudo apt install poppler-utils
fi

# Create absolute path for bin files
#echo "Creating absolute path for tsRNAsearch 'bin' and 'DBs'..."
#myPath=$(pwd)
#sed -i -e "s~ bin~ ${myPath}\/bin~g" bin/tsRNAsearch_single.sh # using tilde as delimiter here instead of slash as myPath variable contains slashes
#sed -i -e "s~ bin~ ${myPath}\/bin~g" tsRNAsearch
#sed -i -e "s~DBs~${myPath}\/DBs~g" bin/tsRNAsearch_single.sh
#sed -i -e "s~DBs~${myPath}\/DBs~g" tsRNAsearch
#sed -i -e "s~bin\/trim~${myPath}\/bin\/trim~g" bin/tsRNAsearch_single.sh
#sed -i -e "s~bin\/feat~${myPath}\/bin\/feat~g" bin/tsRNAsearch_single.sh
#sed -i -e "s~additional~${myPath}\/additional~g" bin/tsRNAsearch_single.sh
#sed -i -e "s~additional~${myPath}\/additional~g" tsRNAsearch

wait # Wait for things to finish running

echo "Done"
