#!/bin/bash

# Check if conda installed:
echo "Checking for conda..."
if [[ $(conda info | wc -l) -ge 1 ]];
then
	echo "Conda is installed"
else
	echo "Conda is not installed. Please install conda (Miniconda2 Linux 64 bit Python2.7) using the following guide: https://docs.conda.io/en/latest/miniconda.html#linux-installer"
	exit 1
fi

### Create conda environment
echo "Creating conda environment for tsRNAearch..."
conda create -y --name tsrnasearch_env python=2.7 # Create new environment with python 2.7
echo -e "

Please activate the environment:

conda activate tsrnasearch_env # Activate new environment

Then run the conda environment installer:

# Choose your species of interest using -s parameter
bash bin/setup_conda_2.sh -s [human/rat/mouse/all]

"
