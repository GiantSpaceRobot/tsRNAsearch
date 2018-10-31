#!/bin/bash

# Usage: ./tiRNA-pipeline.sh -h
# Author: Paul Donovan 
# Email: pauldonovan@rcsi.com
# 19-Oct-2018

asciiArt() { echo "
  __  .____________  _______      _____    __________.__              .__  .__               
_/  |_|__\______   \ \      \    /  _  \   \______   \__|_____   ____ |  | |__| ____   ____  
\   __\  ||       _/ /   |   \  /  /_\  \   |     ___/  \____ \_/ __ \|  | |  |/    \_/ __ \ 
 |  | |  ||    |   \/    |    \/    |    \  |    |   |  |  |_> >  ___/|  |_|  |   |  \  ___/ 
 |__| |__||____|_  /\____|__  /\____|__  /  |____|   |__|   __/ \___  >____/__|___|  /\___  >
                 \/         \/         \/               |__|        \/             \/     \/ 
" 1>&1; }
usage() { echo "Usage (single-end): $0 -o OutputDirectory/ -s SeqFile.fastq.gz
Usage (paired-end): $0 -o OutputDirectory/ -1 SeqFile_1.fastq.gz -2 SeqFile_2.fastq.gz" 1>&2; }
info() { echo "
Options
	
	-h	Print the usage and options information
	-s	Single-end file for analysis. 
	-1	First paired-end file. 
	-2	Second paired-end file.
	-o	Output directory for the results and log files
	
	Input file format should be FASTQ (.fq/.fastq) or gzipped FASTQ (.gz)
	" 1>&2; }

while getopts ":hs:1:2:o:" o; do
    case "${o}" in
		h)
			asciiArt
			usage
			info
			exit
			;;
		s)
			singleFile="$OPTARG"
			;;
		1)
			file1="$OPTARG"
			;;
		2)
			file2="$OPTARG"
			;;
		o)
			outDir="$OPTARG"
			;;
		*)
            echo "Error in input parameters!"
			usage
			exit 1
            ;;
    esac
done

echo "Started at $(date)" # Print pipeline start-time

# Determine if paired-end or single-end files were provided as input
if [ -z "$singleFile" ]; then # If the singleEnd variable is empty
	pairedEnd="True"
else
	pairedEnd="False"
fi

# Check if the output directory exists. If not, create it
if [ ! -d $outDir ]; then 
	mkdir $outDir
	mkdir $outDir/trim_galore_output
	mkdir $outDir/FastQC
fi

# Run Trim_Galore (paired-end or single-end)
if [[ $pairedEnd = "True" ]]; then 
	file1_base=${file1##*/}    # Remove pathname
	basename1="$( cut -d '.' -f 1 <<< "$file1_base" )" # Get filename before first occurence of .	
	suffix1="$( cut -d '.' -f 2- <<< "$file1_base" )" # Get full file suffix/entension
	file2_base=${file2##*/}    # Remove pathname
	basename2="$( cut -d '.' -f 1 <<< "$file2_base" )" # Get filename before first occurence of .	
	suffix2="$( cut -d '.' -f 2- <<< "$file2_base" )" # Get full file suffix/entension
	trim_galore -o $outDir/trim_galore_output/ --paired $file1 $file2  
	printf -v trimmedFile1 "%s_val_1.%s" "$basename1" "$suffix1"
	printf -v trimmedFile2 "%s_val_2.%s" "$basename2" "$suffix2"
	fastqc -o $outDir/FastQC/ -f fastq $outDir/trim_galore_output/$trimmedFile1 $outDir/trim_galore_output/$trimmedFile2 
elif [[ $pairedEnd = "False" ]]; then
	singleFile_base=${singleFile##*/}    # Remove pathname
	singleFile_basename="$( cut -d '.' -f 1 <<< "$singleFile_base" )" # Get filename before first occurence of .	
	suffix="$( cut -d '.' -f 2- <<< "$singleFile_base" )" # Get full file suffix/entension
	trim_galore -o $outDir/FastQC/ $singleFile
	printf -v trimmedFile "%s_trimmed.%s" "$singleFile_basename" "$suffix"
	fastqc -o $outDir/FastQC/ -f fastq $outDir/trim_galore_output/$trimmedFile
fi

# Create a log file with the date, time and name of the input file in it's name
# Run FastQC on files







