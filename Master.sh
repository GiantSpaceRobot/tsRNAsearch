#!/bin/bash
# Author: Paul Donovan 
# Email: pauldonovan@rcsi.com
# 12-12-2018

echo "Started at $(date)"
StartTime="Pipeline initiated at $(date)" 

usage() { echo "
	Usage: $0 -p -d Path/To/Input/Files -o OutputDirectory/ -e ExperimentLayout.csv
	" 1>&2; }
info() { echo "
Options

	-h	Print the usage and options information
	-p	Input files are paired end data
	-d	Directory containing the files for analysis. Directory should have no other contents.
	-o	Output directory for the results and log files
	-e	Optional (but recommended) CSV file containing file names and file groups (see examples in ./additional-files/)
	-t	Number of threads/CPUs to use

	" 1>&2; }

while getopts ":hpt:d:e:o:" o; do
    case "${o}" in
		h)
			asciiArt
			usage
			info
			exit
			;;
		p)
			paired="Yes"
			;;
		d)
			inDir="$OPTARG"
			echo "$inDir"
			;;
		e)
			expFile="$OPTARG"
			echo "$expFile"
			;;
		o)
			outDir="$OPTARG"
			echo "$outDir"
			;;
		t)
			CPUs="$OPTARG"
			echo "$CPUs"
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
    echo "No command line parameters provided!"
	usage
    exit 0
fi

for f in $inDir/*; do
	mkdir -p $outDir
	mkdir -p $outDir/Results
	mkdir -p $outDir/Results/Data
	mkdir -p $outDir/Results/Plots
	file_base=$(basename $f)
	#filename=${file_base%.*}
	#singleFile_base=${singleFile##*/}    # Remove pathname
	filename="$( cut -d '.' -f 1 <<< "$file_base" )" 
	./tiRNA-pipeline.sh -s "$f" -o "$outDir/Results/$filename" -p "$CPUs" -T
	cp $outDir/Results/$filename/Plots/* $outDir/Results/Plots/
	wait
	cat $outDir/Results/$filename/HTSeq-count-output/*.count | grep -v ^__ | sort -k1,1 > $outDir/Results/Data/$filename.all_features.count
	sed -i '1s/^/Features\t'"$filename"'\n/' $outDir/Results/Data/$filename.all_features.count # Add column headers
	readsMapped=$(awk '{sum+=$2} END{print sum;}' $outDir/Results/Data/$filename.all_features.count)
	echo "Reads mapped: $readsMapped" >> $outDir/Results/$filename/Stats.log
	#if [ ! "$expFile" ];
	#	echo -e "$filename" >> $outDir/Results/Plots/FilenamesForR.txt
	#
	#fi
done

### Gather count files
awk '{print $1}' $outDir/Results/Data/$filename.all_features.count > $outDir/Results/Data/HTSeq.all_features
for f in $outDir/Results/Data/*count; do
	awk '{print $2}' $f | paste $outDir/Results/Data/HTSeq.all_features - >> $outDir/Results/Data/HTSeq.temp
	mv $outDir/Results/Data/HTSeq.temp $outDir/Results/Data/HTSeq.all_features
done

mv $outDir/Results/Data/HTSeq.all_features $outDir/Results/Plots/HTSeq.all_features.count

### Determine if experiment layout file was provided or not. If not, ry and figure out which files group together using R.
myPath=$(pwd) #Get working dir to recreate full path for R script execution
if [ ! "$expFile" ]; then
    echo "No experiment layout plan provided. This will now be created prior to the formal DESeq2 analysis."
    Rscript --vanilla scripts/DESeq2_tiRNA-pipeline.R "$myPath/$outDir/Results/Data/"
else
    echo "An experiment layout plan was provided. Carrying out DESeq2 analysis now."
    Rscript --vanilla scripts/DESeq2_tiRNA-pipeline.R "$expFile" "$myPath/$outDir/Results/Data/"
fi






#rm -rf Data

echo "Finished at $(date)"

