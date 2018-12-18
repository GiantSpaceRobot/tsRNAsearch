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
	mkdir -p $outDir/Results/temp_files
	file_base=$(basename $f)
	filename=${file_base%.*}
	./tiRNA-pipeline.sh -s "$f" -o "$outDir/Results/$filename" -p "$CPUs" -T
	wait
	cat $outDir/Results/$filename/HTSeq-count-output/*.count | grep -v ^__ | sort -k1,1 > $outDir/Results/temp_files/$filename.all_features.count
	sed -i '1s/^/Features\t'"$filename"'\n/' $outDir/Results/temp_files/$filename.all_features.count
	readsMapped=$(awk '{sum+=$2} END{print sum;}' $outDir/Results/temp_files/$filename.all_features.count)
	echo $readsMapped >> $outDir/Results/$filename/Stats.log
done

### Gather count files
awk '{print $1}' $outDir/Results/temp_files/$filename.all_features.count > $outDir/Results/temp_files/HTSeq.all_features
for f in $outDir/Results/temp_files/*count; do
	awk '{print $2}' $f | paste $outDir/Results/temp_files/HTSeq.all_features - >> $outDir/Results/temp_files/HTSeq.temp
	mv $outDir/Results/temp_files/HTSeq.temp $outDir/Results/temp_files/HTSeq.all_features
done

mv $outDir/Results/temp_files/HTSeq.all_features $outDir/Results/HTSeq.all_features.count


### Determine if experiment layout file was provided or not. Try and figure out which files group together using R. 
if [ ! "$expFile" ]; then
	Rscript --vanilla DESeq2_tiRNA-pipeline.R $inDir/
else
	Rscript --vanilla DESeq2_tiRNA-pipeline.R $expFile $inDir/
fi





#rm -rf temp_files

echo "Finished at $(date)"

