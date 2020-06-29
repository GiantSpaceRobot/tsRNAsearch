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
	STAR --runThreadN $threads --runMode genomeGenerate --genomeDir DBs/species_index/human-ncRNAs/ --genomeFastaFiles DBs/human_tRNAs-and-ncRNAs_relative_cdhit.fa --genomeSAindexNbases 8
}
function mouse_setup () {
	mkdir -p DBs/species_index/mouse-ncRNAs
	STAR --runThreadN $threads --runMode genomeGenerate --genomeDir DBs/species_index/mouse-ncRNAs/ --genomeFastaFiles DBs/mouse_tRNAs-and-ncRNAs_relative_cdhit.fa
}
function rat_setup () {
	mkdir -p DBs/species_index/rat-ncRNAs
	STAR --runThreadN $threads --runMode genomeGenerate --genomeDir DBs/species_index/rat-ncRNAs/ --genomeFastaFiles DBs/rat_tRNAs-and-ncRNAs_relative_cdhit.fa
}

### Set up conda environment
conda info # Check if conda installed
conda create -y --name tsrnasearch_env python=2.7 # Create new environment with python 2.7
source activate tsrnasearch_env # Activate new environment

### install all required tools and packages
conda install -y -c bioconda star
conda install -y -c bioconda trim-galore
conda install -y numpy
conda install -y -c r r # Install R
conda install -y -c r r-essentials
conda install -y -c conda-forge r-metap
conda install -y -c bioconda bioconductor-deseq2
conda install -y -c conda-forge r-ggrepel
conda install -y -c conda-forge r-gplots
conda install -y -c conda-forge r-venndiagram
conda install -y -c bioconda bioconductor-genomeinfodb
conda install -y -c bioconda bioconductor-enhancedvolcano
conda install -y -c bioconda samtools

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

source conda deactivate

echo "Run 'source activate '"
