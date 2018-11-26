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
	-s	Single-end file for analysis 
	-1	First paired-end file
	-2	Second paired-end file
	-o	Output directory for the results and log files
	-T	Use Tophat2 to carry out alignment steps (default: HISAT2)

	Input file format should be FASTQ (.fq/.fastq) or gzipped FASTQ (.gz)
	" 1>&2; }

while getopts ":hTs:1:2:o:" o; do
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
		T)
			oldAligner="yes"
			;;
		*)
            echo "Error in input parameters!"
			usage
			exit 1
            ;;
    esac
done

echo "Started at $(date)" # Print pipeline start-time

# If I want to use Naser's tRNA DB, I can change the value of a variable to mine or his without adding more code 

# Function to pad text with characters to make sentences stand out more
function string_padder () {
	string=$1
	count=${#string}
	v=$(printf "%-${count}s" "=")
	padding="${v// /=}"
	echo "

	$padding
	$string
	$padding

	"
}

# Determine number of CPUs available
maxCPUs=$(grep -c ^processor /proc/cpuinfo)

# If the variable "CPUs" is not an integer, default to use a value of 1 
re='^[0-9]+$'
if ! [[ $maxCPUs =~ $re ]] ; then
	echo "Error: Variable 'maxCPUs' is not an integer" >&2
	maxCPUs=1
fi
CPUs=$(echo "$maxCPUs * 0.75" | bc | awk '{print int($1+0.5)}') # Use 75% of max CPU number

# Determine if paired-end or single-end files were provided as input
if [ -z "$singleFile" ]; then # If the singleEnd variable is empty
	pairedEnd="True"
else
	pairedEnd="False"
fi

# Check if the output directory exists. If not, create it
message="Creating directory structure"
string_padder "$message"
if [ ! -d $outDir ]; then 
	mkdir $outDir
	mkdir $outDir/trim_galore_output
	mkdir $outDir/FastQC
	mkdir $outDir/tRNA-alignment
	mkdir $outdir/snomiRNA-alignment
fi

# Run Trim_Galore (paired-end or single-end)
message="Trimming reads using Trim Galore"
string_padder "$message"
if [[ $pairedEnd = "True" ]]; then 

	echo "
	Input: paired-end read files
	"

	file1_base=${file1##*/}    # Remove pathname
	basename1="$( cut -d '.' -f 1 <<< "$file1_base" )" # Get filename before first occurence of .	
	suffix1="$( cut -d '.' -f 2- <<< "$file1_base" )" # Get full file suffix/entension
	file2_base=${file2##*/}    # Remove pathname
	basename2="$( cut -d '.' -f 1 <<< "$file2_base" )" # Get filename before first occurence of .	
	suffix2="$( cut -d '.' -f 2- <<< "$file2_base" )" # Get full file suffix/entension
	
	# Run Trim_Galore on paired-end read files
	trim_galore --stringency 10 --length 15 -o $outDir/trim_galore_output/ --paired $file1 $file2  
	printf -v trimmedFile1 "%s_val_1.%s" "$basename1" "$suffix1"
	printf -v trimmedFile2 "%s_val_2.%s" "$basename2" "$suffix2"
	
	message="Trimming complete. Moving on to FastQC analysis"
	string_padder "$message"

	# Run FastQC on newly trimmed files
	fastqc -o $outDir/FastQC/ -f fastq $outDir/trim_galore_output/$trimmedFile1 $outDir/trim_galore_output/$trimmedFile2 
	message="Finished running FastQC. Moving onto alignment to tRNA database"
	string_padder "$message"
	
	# Align trimmed reads to tRNA database using HISAT2/Tophat2
	if [[ $oldAligner = "yes" ]]; then
		tophat2 -p $CPUs -x 1 -o $outDir/tRNA-alignment/ DBs/bowtie2_index/hg19-wholetRNA-CCA $outDir/trim_galore_output/$trimmedFile1 $outDir/trim_galore_output/$trimmedFile2
		# Tophat2 using HG38 genome below:
		#tophat2 -p $CPUs -x 1 -o $outDir/tRNA-alignment/ DBs/bowtie2_index/hg38-tRNAs_CCA $outDir/trim_galore_output/$trimmedFile1 $outDir/trim_galore_output/$trimmedFile2 
		bedtools bamtofastq -i $outDir/tRNA-alignment/unmapped.bam -fq $outDir/tRNA-alignment/unmapped_1.fastq -fq2 $outDir/tRNA-alignment/unmapped_2.fastq
	else
		echo $trimmedFile1
		hisat2 -p $CPUs -x DBs/hisat2_index/hg38-tRNAs_CCA -1 $outDir/trim_galore_output/$trimmedFile1 -2 $outDir/trim_galore_output/$trimmedFile2 -S $outDir/tRNA-alignment/aligned_tRNAdb.sam
	fi

elif [[ $pairedEnd = "False" ]]; then

	echo "
	Input: single-end read file
	"
	singleFile_base=${singleFile##*/}    # Remove pathname
	singleFile_basename="$( cut -d '.' -f 1 <<< "$singleFile_base" )" # Get filename before first occurence of .	
	suffix="$( cut -d '.' -f 2- <<< "$singleFile_base" )" # Get full file suffix/entension

	# Run Trim_Galore on single read file
	trim_galore --stringency 10 --length 15 -o $outDir/trim_galore_output/ $singleFile
	printf -v trimmedFile "%s_trimmed.%s" "$singleFile_basename" "$suffix"
	message="Trimming complete. Moving on to FastQC analysis"
	string_padder "$message"

	# Run FastQC on newly trimmed file
	fastqc -o $outDir/FastQC/ -f fastq $outDir/trim_galore_output/$trimmedFile
	message="Finished running FastQC. Moving onto alignment to tRNA database"
	string_padder "$message"

	# Align trimmed reads to tRNA database using HISAT2/Tophat2
	if [[ $oldAligner = "yes" ]]; then
		tophat2 -p $CPUs -x 1 -o $outDir/tRNA-alignment/ DBs/bowtie2_index/hg19-wholetRNA-CCA $outDir/trim_galore_output/$trimmedFile
		# Tophat2 using HG38 genome below
		#tophat2 -p $CPUs -x 1 -o $outDir/tRNA-alignment/ DBs/bowtie2_index/hg38-tRNAs_CCA $outDir/trim_galore_output/$trimmedFile
		bedtools bamtofastq -i $outDir/tRNA-alignment/unmapped.bam -fq $outDir/tRNA-alignment/unmapped.fastq
		tophat2 -p $CPUs -x 1 -o $outDir/snomiRNA-alignment/ DBs/bowtie2_index/hg19-snomiRNA $outDir/tRNA-alignment/unmapped.fastq
		# tophat2 has a bug that causes it to crash if none of the reads map (I think). If it crashes (and no results are generated), use the unmapped.fastq from the tRNA alignment step
		if [ ! -f $outDir/snomiRNA-alignment/accepted_hits.bam ]; then
			echo "snomiRNA accepted_hits.bam file not found! Using unmapped.fastq file from tRNA alignment output"
			unmappedReads="$outDir/tRNA-alignment/unmapped.fastq"
		else
			bedtools bamtofastq -i $outDir/snomiRNA-alignment/unmapped.bam -fq $outDir/snomiRNA-alignment/unmapped.fastq
			unmappedReads="$outDir/snomiRNA-alignment/unmapped.fastq"
		fi
		tophat2 -p $CPUs -o $outDir/ncRNA-mRNA-alignment/ DBs/bowtie2_index/### $unmappedReads

	else
		hisat2 -p $CPUs -x DBs/hisat2_index/hg38-tRNAs_CCA -U $outDir/trim_galore_output/$trimmedFile -S $outDir/tRNA-alignment/aligned_tRNAdb.sam
	fi
fi

# Create a log file with the date, time and name of the input file in it's name
# Run FastQC on files







