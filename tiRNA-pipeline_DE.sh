#!/bin/bash
# Author: Paul Donovan 
# Email: pauldonovan@rcsi.com
# 12-12-2018

asciiArt() { echo "

	  __  .____________  _______      _____    __________.__              .__  .__               
	_/  |_|__\______   \ \      \    /  _  \   \______   \__|_____   ____ |  | |__| ____   ____  
	\   __\  ||       _/ /   |   \  /  /_\  \   |     ___/  \____ \_/ __ \|  | |  |/    \_/ __ \ 
	 |  | |  ||    |   \/    |    \/    |    \  |    |   |  |  |_> >  ___/|  |_|  |   |  \  ___/ 
	 |__| |__||____|_  /\____|__  /\____|__  /  |____|   |__|   __/ \___  >____/__|___|  /\___  >
			 \/         \/         \/               |__|        \/             \/     \/ 
	
	" 1>&1; }
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
	-t	Number of threads/CPUs to use {default is to calculate the number of processors and use 75%}

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
			info
			exit 1
            ;;
    esac
done

### If no command line arguments provided, quit
if [ -z "$*" ] ; then
    echo "No command line parameters provided!"
	asciiArt
	usage
    info
	exit 1
fi

### If the pathname specified by $inDir does not begin with a slash, quit (we need full path name)
if [[ ! $inDir = /* ]]; then
	echo "Error: File paths must absolute. Please specify the full path for the input directory."
	exit 1
fi

### If the pathname specified by $outDir does not begin with a slash, quit (we need full path name)
if [[ ! $outDir = /* ]]; then
	echo "Error: File paths must absolute. Please specify the full path for the output directory."
	exit 1
fi

### If the pathname specified by $expFile does not begin with a slash, quit (we need full path name)
if [ "$expFile" ]; then
    if [[ ! $expFile = /* ]]; then
        echo "Error: File paths must absolute. Please specify the full path for the experiment layout file."
        exit 1
    fi
fi

echo "Started at $(date)"
StartTime="Pipeline initiated at $(date)"

for f in $inDir/*; do
	mkdir -p $outDir
	mkdir -p $outDir/Results
	mkdir -p $outDir/Results/Data
	mkdir -p $outDir/Results/Plots
	file_base=$(basename $f)
	filename="$( cut -d '.' -f 1 <<< "$file_base" )" 
	./tiRNA-pipeline.sh -s "$f" -o "$outDir/Results/$filename" -p "$CPUs" -T
	cp $outDir/Results/$filename/Data_and_Plots/* $outDir/Results/Plots/
	wait
	cat $outDir/Results/$filename/HTSeq-count-output/*.count | grep -v ^__ | sort -k1,1 > $outDir/Results/Data/$filename.all_features.count
	sed -i '1s/^/Features\t'"$filename"'\n/' $outDir/Results/Data/$filename.all_features.count # Add column headers
	readsMapped=$(awk '{sum+=$2} END{print sum;}' $outDir/Results/Data/$filename.all_features.count)
	echo "Reads mapped: $readsMapped" >> $outDir/Results/$filename/Stats.log
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
    condition1=$( awk -F ',' 'NR == 1 {print $2}' $myPath/$outDir/Results/Data/predicted_exp_layout.csv ) # Get element in first row second column (condition)
    grep $condition1 $myPath/$outDir/Results/Data/predicted_exp_layout.csv > $myPath/$outDir/Results/Data/predicted_exp_layout_cond1.csv
    condition2=$( grep -v $condition1 $myPath/$outDir/Results/Data/predicted_exp_layout.csv | awk -F ',' 'NR == 1 {print $2}') # Get the second condition using the inverse of the first one
    grep $condition2 $myPath/$outDir/Results/Data/predicted_exp_layout.csv > $myPath/$outDir/Results/Data/predicted_exp_layout_cond2.csv
    ### Concatenate the genomecov files for the replicates of condition 1
    count=0
    for fname in $(cat $myPath/$outDir/Results/Data/predicted_exp_layout_cond1.csv); do
        count=$((count + 1))
        cond1="$(cut -d ',' -f1 <<< $fname)"
        cond1base="$(cut -d '.' -f1 <<< $cond1)"
        find $myPath/$outDir/ -type f -name "$cond1base*tiRNA.genomecov" -exec cp {} $myPath/$outDir/Results/Data/condition1_file$count.genomecov \;
    done
    paste $myPath/$outDir/Results/Data/condition1_file*.genomecov > $myPath/$outDir/Results/Data/condition1_concatenated.genomecov
    ### Concatenate the genomecov files for the replicates of condition 2 
    count=0
    for fname in $(cat $myPath/$outDir/Results/Data/predicted_exp_layout_cond2.csv); do
        count=$((count + 1))
        cond2="$(cut -d ',' -f1 <<< $fname)"
        cond2base="$(cut -d '.' -f1 <<< $cond2)"
        find $myPath/$outDir/ -type f -name "$cond2base*tiRNA.genomecov" -exec cp {} $myPath/$outDir/Results/Data/condition2_file$count.genomecov \;
    done
    paste $myPath/$outDir/Results/Data/condition2_file*.genomecov > $myPath/$outDir/Results/Data/condition2_concatenated.genomecov
    ### Get mean and standard deviation for both files
    Rscript scripts/Mean_Stdev.R $myPath/$outDir/Results/Data/condition1_concatenated.genomecov $myPath/$outDir/Results/Data/condition1_concatenated_mean_stdev.genomecov
    Rscript scripts/Mean_Stdev.R $myPath/$outDir/Results/Data/condition2_concatenated.genomecov $myPath/$outDir/Results/Data/condition2_concatenated_mean_stdev.genomecov
    ### Plot DEGs arg1 and 2 are inputs, arg 3 is output pdf, arg 4 is mean coverage cutoff (plot features with coverage above this)
    Rscript scripts/Bedgraph_plotter_DEGs.R $myPath/$outDir/Results/Data/condition1_concatenated_mean_stdev.genomecov $myPath/$outDir/Results/Data/condition2_concatenated_mean_stdev.genomecov $myPath/$outDir/Results/Plots/${condition1}_vs_${condition2}_DEGs.pdf 20

else
    echo "An experiment layout plan was provided. Carrying out DESeq2 analysis now."
    Rscript --vanilla scripts/DESeq2_tiRNA-pipeline.R "$expFile" "$myPath/$outDir/Results/Data/"
    condition1=$( awk -F ',' 'NR == 1 {print $2}' "$expFile" ) # Get element in first row second column (condition)
    grep $condition1 "$expFile" > $myPath/$outDir/Results/Data/predicted_exp_layout_cond1.csv
    condition2=$( grep -v $condition1 "$expFile" | awk -F ',' 'NR == 1 {print $2}') # Get the second condition using the inverse of the first one
    grep $condition2 "$expFile" > $myPath/$outDir/Results/Data/predicted_exp_layout_cond2.csv
    ### Concatenate the genomecov files for the replicates of condition 1
    count=0
    for fname in $(cat $myPath/$outDir/Results/Data/predicted_exp_layout_cond1.csv); do
        count=$((count + 1))
        cond1="$(cut -d ',' -f1 <<< $fname)"
        cond1base="$(cut -d '.' -f1 <<< $cond1)"
        find $myPath/$outDir/ -type f -name "$cond1base*tiRNA.genomecov" -exec cp {} $myPath/$outDir/Results/Data/condition1_file$count.genomecov \;
    done
    paste $myPath/$outDir/Results/Data/condition1_file*.genomecov > $myPath/$outDir/Results/Data/condition1_concatenated.genomecov
    ### Concatenate the genomecov files for the replicates of condition 2 
    count=0
    for fname in $(cat $myPath/$outDir/Results/Data/predicted_exp_layout_cond2.csv); do
        count=$((count + 1))
        cond2="$(cut -d ',' -f1 <<< $fname)"
        cond2base="$(cut -d '.' -f1 <<< $cond2)"
        find $myPath/$outDir/ -type f -name "$cond2base*tiRNA.genomecov" -exec cp {} $myPath/$outDir/Results/Data/condition2_file$count.genomecov \;
    done
    paste $myPath/$outDir/Results/Data/condition2_file*.genomecov > $myPath/$outDir/Results/Data/condition2_concatenated.genomecov
    ### Get mean and standard deviation for both files
    Rscript scripts/Mean_Stdev.R $myPath/$outDir/Results/Data/condition1_concatenated.genomecov $myPath/$outDir/Results/Data/condition1_concatenated_mean_stdev.genomecov
    Rscript scripts/Mean_Stdev.R $myPath/$outDir/Results/Data/condition2_concatenated.genomecov $myPath/$outDir/Results/Data/condition2_concatenated_mean_stdev.genomecov
    ### Plot DEGs arg1 and 2 are inputs, arg 3 is output pdf, arg 4 is mean coverage cutoff (plot features with coverage above this)
    Rscript scripts/Bedgraph_plotter_DEGs.R $myPath/$outDir/Results/Data/condition1_concatenated_mean_stdev.genomecov $myPath/$outDir/Results/Data/condition2_concatenated_mean_stdev.genomecov $myPath/$outDir/Results/Plots/${condition1}_vs_${condition2}_DEGs.pdf 20
fi


#rm -rf Data

echo "Finished at $(date)"

