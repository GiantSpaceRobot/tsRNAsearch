#!/bin/bash
# Author: Paul Donovan 
# Email: pauldonovan@rcsi.com
# 12-12-2018

asciiArt() { echo '
  
 _        _____ _   _   ___                          _     
| |      | ___ \ \ | | / _ \                        | |    
| |_ ___ | |_/ /  \| |/ /_\ \___  ___  __ _ _ __ ___| |__  
| __/ __||    /| . ` ||  _  / __|/ _ \/ _` | `__/ __| `_ \ 
| |_\__ \| |\ \| |\  || | | \__ \  __/ (_| | | | (__| | | |
 \__|___/\_| \_\_| \_/\_| |_/___/\___|\__,_|_|  \___|_| |_|
                                                          
   
   ' 1>&1; }
usage() { echo "
	Usage: $0 -g human/mouse -d Path/To/Input/Files -o OutputDirectory/ -e ExperimentLayout.csv -t CPU-number
	" 1>&2; }
info() { echo "
Options

	-h	Print the usage and options information
	-g	Analyse datasets against 'human' or 'mouse'? {default: human}
	-d	Directory containing the files for analysis. Directory should have no other contents.
	-o	Output directory for the results and log files
	-e	CSV file containing file names and file groups (see examples in additional-files/)
	-t	Number of threads to use {default is to calculate the number of processors and use 75%}
	-A	Plot all features? yes/no {default: no (only plot differentially expressed features)}

	" 1>&2; }

while getopts ":hg:A:t:d:e:o:" o; do
    case "${o}" in
		h)
			asciiArt
			usage
			info
			exit
			;;
		g)
			genome="$OPTARG"
			;;
		d)
			inDir="$OPTARG"
			;;
		e)
			expFile="$OPTARG"
			;;
		o)
			outDir="$OPTARG"
			;;
		t)
			threads="$OPTARG"
			;;
		A)
			Plots="$OPTARG"
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
#if [[ ! $inDir = /* ]]; then
#	echo "Error: File paths must absolute. Please specify the full path for the input directory."
#	exit 1
#fi

### If the pathname specified by $outDir does not begin with a slash, quit (we need full path name)
#if [[ ! $outDir = /* ]]; then
#	echo "Error: File paths must absolute. Please specify the full path for the output directory."
#	exit 1
#fi

### If the pathname specified by $expFile does not begin with a slash, quit (we need full path name)
#if [ "$expFile" ]; then
#    if [[ ! $expFile = /* ]]; then
#        echo "Error: File path must absolute. Please specify the full path for the experiment layout file."
#        exit 1
#    fi
#fi

### Get working dir to recreate full path for R script execution
myPath=$(pwd) 

### Define functions
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

function DataTransformations () {
### Given half of a layout file, transform data
count=0
for fname in $(cat $1); do
	count=$((count + 1))
	cond1="$(cut -d ',' -f1 <<< $fname)"
	cond1base="$(cut -d '.' -f1 <<< $cond1)"
	echo "Finding file for $fname ($cond1base)..."
	find $myPath/$outDir/ -type f -name "$cond1base*tsRNA.depth" -exec cp {} $myPath/$outDir/Data/Intermediate-files/$2_file$count.tsRNA.depth \; & # Gather tsRNAs
	find $myPath/$outDir/ -type f -name "$cond1base*snomiRNA.depth" -exec cp {} $myPath/$outDir/Data/Intermediate-files/$2_file$count.snomiRNA.depth \; & # Gather sno/miRNAs
	wait
	mapped=$(grep "mapped" $myPath/$outDir/$cond1base/Stats.log | awk '{print $3}')
	echo "Converting raw counts to RPM..."
	python bin/Depth-to-Depth_RPM.py "$myPath/$outDir/Data/Intermediate-files/$2_file$count.tsRNA.depth" "$mapped" "$myPath/$outDir/Data/Intermediate-files/$2_file$count.tsRNA.depth.readspermil" &
	python bin/Depth-to-Depth_RPM.py "$myPath/$outDir/Data/Intermediate-files/$2_file$count.snomiRNA.depth" "$mapped" "$myPath/$outDir/Data/Intermediate-files/$2_file$count.snomiRNA.depth.readspermil" &
	### Left-over reads:
	#cp $myPath/$outDir/$cond1base/tRNA-alignment/$cond1base.tRNAs-almost-mapped_RPM.depth $myPath/$outDir/Data/Intermediate-files/$2.$cond1base.tRNAs-almost-mapped_RPM.depth
done
wait
### Concatenate data horizontally
paste $myPath/$outDir/Data/Intermediate-files/$2_file*.tsRNA.depth.readspermil > $myPath/$outDir/Data/Intermediate-files/tsRNA.$2_concatenated.depth &
paste $myPath/$outDir/Data/Intermediate-files/$2_file*.snomiRNA.depth.readspermil > $myPath/$outDir/Data/Intermediate-files/snomiRNA.$2_concatenated.depth &
#paste $myPath/$outDir/Data/Intermediate-files/$2.*.tRNAs-almost-mapped_RPM.depth > $myPath/$outDir/Data/Intermediate-files/Combined.$2.tRNAs-almost-mapped_RPM.depth & # leftover reads
### Concatenate data vertically
cat $myPath/$outDir/Data/Intermediate-files/$2_file*.tsRNA.depth.readspermil | sort -k1,1 -k2,2n > $myPath/$outDir/Data/Intermediate-files/tsRNA.$2_concatenated.depthVert &
cat $myPath/$outDir/Data/Intermediate-files/$2_file*.snomiRNA.depth.readspermil | sort -k1,1 -k2,2n > $myPath/$outDir/Data/Intermediate-files/snomiRNA.$2_concatenated.depthVert &
wait
### Leftovers:
#python bin/MeanCalculator.py $myPath/$outDir/Data/Intermediate-files/$myPath/$outDir/Data/Intermediate-files/Combined.$2.tRNAs-almost-mapped_RPM.depth $myPath/$outDir/Data/Intermediate-files/DataTransformations/Combined.$2.tRNAs-almost-mapped_RPM.depth.mean &
#sort -k1,1 -k2,2n $myPath/$outDir/Data/Intermediate-files/DataTransformations/Combined.$2.tRNAs-almost-mapped_RPM.depth.mean > $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_Combined.$2.tRNAs-almost-mapped_RPM.depth.mean &
### Calculate mean
python bin/MeanCalculator.py $myPath/$outDir/Data/Intermediate-files/tsRNA.$2_concatenated.depth $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.$2_concatenated.depth.mean &
python bin/MeanCalculator.py $myPath/$outDir/Data/Intermediate-files/snomiRNA.$2_concatenated.depth $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.$2_concatenated.depth.mean &
wait
### Sort output
sort -k1,1 -k2,2n $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.$2_concatenated.depth.mean > $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_tsRNA.$2_concatenated.depth.mean &
sort -k1,1 -k2,2n $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.$2_concatenated.depth.mean > $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_snomiRNA.$2_concatenated.depth.mean &
wait
}

                 #############################
                 ##### Start of pipeline #####
                 #############################

pipeline_start="Started project analysis on $date"
string_padder "$pipeline_start"
StartTime="Pipeline initiated at $(date)"

### Are we analysing Human or Mouse? -g parameter
if [ "$genome" ]; then
	if [[ $genome == "human" ]]; then
		snomiRNAGTF="DBs/hg19-snomiRNA_cdhit.gtf"
	elif [[ $genome == "mouse" ]]; then
		snomiRNAGTF="DBs/mouse_snomiRNAs_relative_cdhit.gtf"
	fi
else
	snomiRNAGTF="DBs/hg19-snomiRNA_cdhit.gtf"
	genome="human"
fi

### Print parameters used
echo -e "Parameters:
	Genome (-g): $genome
	Input directory containing fastq/fastq.gz files (-d): $indir
	Experiment layout file (-e): $expFile
	Output directory that tsRNAsearch will create and populate with results (-o): $outDir
	Number of threads to use for the analysis (-t): $threads
	"

### If -A parameter was not provided, default is to only plot differentially expressed features
if [ ! "$Plots" ]; then
	Plots="no"
else
	Plots="yes"
fi

### Create dir substructure
for f in $inDir/*; do
	mkdir -p $outDir
	mkdir -p $outDir/Data
	mkdir -p $outDir/Data/Intermediate-files
	mkdir -p $outDir/Data/Intermediate-files/DataTransformations
	mkdir -p $outDir/Plots
	file_base=$(basename $f)
	filename="$( cut -d '.' -f 1 <<< "$file_base" )" 
	analysis="Beginning analysis of $filename using tsRNAsearch"
	string_padder $analysis
	tsRNAsearch.sh -g "$genome" -s "$f" -o "$outDir/$filename" -t "$threads" -A "$Plots"
	wait
	cat $outDir/$filename/FCount-count-output/*.count | grep -v ^__ | sort -k1,1 > $outDir/Data/Intermediate-files/$filename.all-features.count	
	readsMapped=$(awk '{sum+=$2} END{print sum;}' $outDir/Data/Intermediate-files/$filename.all-features.count)
	cp $outDir/$filename/Data_and_Plots/$filename.all-features.rpm.count $outDir/Data/	
	cp $outDir/$filename/FCount-count-output/$filename.all-features.count $outDir/Data/Intermediate-files/ 
done

### Gather raw count files
awk '{print $1}' $outDir/Data/Intermediate-files/$filename.all-features.count > $outDir/Data/Intermediate-files/FCount.all-features # Get feature names
for f in $outDir/Data/Intermediate-files/*count; do
	awk '{print $2}' $f | paste $outDir/Data/Intermediate-files/FCount.all-features - >> $outDir/Data/Intermediate-files/FCount.temp
	mv $outDir/Data/Intermediate-files/FCount.temp $outDir/Data/Intermediate-files/FCount.all-features
done
mv $outDir/Data/Intermediate-files/FCount.all-features $outDir/Data/FCount.all-features.raw.count

### Carry out DESeq2 analysis
string_padder "Carrying out DESeq2 analysis..."
condition1=$( awk -F ',' 'NR == 1 {print $2}' "$expFile" ) # Get element in first row second column (condition)
condition2=$( grep -v $condition1 "$expFile" | awk -F ',' 'NR == 1 {print $2}') # Get the second condition using the inverse of the first one
Rscript --vanilla bin/DESeq2_tsRNAsearch.R "$myPath/$outDir/Data/Intermediate-files/" "${condition1}_vs_${condition2}" "$snomiRNAGTF" "$expFile"
grep $condition1 "$expFile" > $myPath/$outDir/Data/Intermediate-files/predicted_exp_layout_cond1.csv
grep $condition2 "$expFile" > $myPath/$outDir/Data/Intermediate-files/predicted_exp_layout_cond2.csv

### Transform data 
DataTransformations $myPath/$outDir/Data/Intermediate-files/predicted_exp_layout_cond1.csv ${condition1} 
DataTransformations $myPath/$outDir/Data/Intermediate-files/predicted_exp_layout_cond2.csv ${condition2}

### Get distribution scores (standard deviation of RPM difference between samples multiplied 
### by standard deviation of percent difference between samples, divided by 1000) of the features
string_padder "tsRNAs: Generating Distribution Scores..."

### Concatenate horizontally to make a dataframe
paste $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_tsRNA.${condition1}_concatenated.depth.mean $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_tsRNA.${condition2}_concatenated.depth.mean > $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.cond1-vs-cond2.mean
### Calculate StdDev
python bin/Mean-to-RelativeDifference.py $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.cond1-vs-cond2.mean $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.cond1-vs-cond2.stddev
### Calculate distribution scores
Rscript bin/DistributionScore.R $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.cond1-vs-cond2.stddev $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.cond1-vs-cond2 $snomiRNAGTF
### Sort the output but not the header
cat $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.cond1-vs-cond2.high-distribution-score.txt | awk 'NR<2{print $0;next}{print $0| "sort -k5,5nr"}' > $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.cond1-vs-cond2.high-distribution-score.sorted.txt
cat $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.cond1-vs-cond2.all-features.txt | awk 'NR<2{print $0;next}{print $0| "sort -k11,11nr"}' > $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.cond1-vs-cond2.all-features.sorted.txt
### Repeat these steps for reads mapped to tRNA groups (multi-mapping reads)
#paste $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_Combined.${condition1}.tRNAs-almost-mapped_RPM.depth.mean $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_Combined.${condition2}.tRNAs-almost-mapped_RPM.depth.mean > $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.multi-mappers.cond1-vs-cond2.mean
#python bin/Mean-to-RelativeDifference.py $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.multi-mappers.cond1-vs-cond2.mean $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.multi-mappers.cond1-vs-cond2.stddev
#Rscript bin/DistributionScore.R $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.multi-mappers.cond1-vs-cond2.stddev $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.multi-mappers.cond1-vs-cond2 $snomiRNAGTF
#cat $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.multi-mappers.cond1-vs-cond2.high-distribution-score.txt | awk 'NR<2{print $0;next}{print $0| "sort -k5,5nr"}' > $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.multi-mappers.cond1-vs-cond2.high-distribution-score.sorted.txt
#cat $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.multi-mappers.cond1-vs-cond2.all-features.txt | awk 'NR<2{print $0;next}{print $0| "sort -k11,11nr"}' > $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.multi-mappers.cond1-vs-cond2.all-features.sorted.txt
###

string_padder "snoRNAs/miRNAs: Generating Distribution Scores..."
### snomiRNAs
### Concatenate horizontally to make a dataframe
paste $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_snomiRNA.${condition1}_concatenated.depth.mean $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_snomiRNA.${condition2}_concatenated.depth.mean > $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.mean
### Calculate StdDev
python bin/Mean-to-RelativeDifference.py $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.mean $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.stddev
### Calculate distribution score
Rscript bin/DistributionScore.R $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.stddev $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2 $snomiRNAGTF
### Sort the output but not the header
cat $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.high-distribution-score.txt | awk 'NR<2{print $0;next}{print $0| "sort -k6,6nr"}' > $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.high-distribution-score.sorted.txt
### Copy txt files illustrating features with very different distributions to the Results directory:
cp $myPath/$outDir/Data/Intermediate-files/DataTransformations/tsRNA.cond1-vs-cond2.high-distribution-score.sorted.txt $myPath/$outDir/Data/${condition1}_vs_${condition2}_High-distribution-tsRNAs.txt 
cp $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.high-distribution-score.sorted.txt $myPath/$outDir/Data/${condition1}_vs_${condition2}_High-distribution-snomiRNAs.txt


string_padder "Generating Cleavage Scores..." 
### Test whether certain features are cleaved in one condition vs the other
Rscript bin/CleavageScore.R $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_tsRNA.${condition1}_concatenated.depth.mean $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_tsRNA.${condition2}_concatenated.depth.mean $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_tsRNAs $snomiRNAGTF &
Rscript bin/CleavageScore.R $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_snomiRNA.${condition1}_concatenated.depth.mean $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_snomiRNA.${condition2}_concatenated.depth.mean $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_snomiRNAs $snomiRNAGTF &
wait
### Sort the output but not the header
cat $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_tsRNAs.high-cleavage-score.txt | awk 'NR<2{print $0;next}{print $0| "sort -k11,11nr"}' > $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_High-cleavage-tsRNAs.txt
cat $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_snomiRNAs.high-cleavage-score.txt | awk 'NR<2{print $0;next}{print $0| "sort -k11,11nr"}' > $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_High-cleavage-snomiRNAs.txt
### Copy sorted high distribution files to the Data dir
cp $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_High-cleavage-snomiRNAs.txt $myPath/$outDir/Data/
cp $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_High-cleavage-tsRNAs.txt $myPath/$outDir/Data/
### Copy PDF files to Plots dir
cp $myPath/$outDir/Data/Intermediate-files/*.pdf $myPath/$outDir/Plots/
cp $myPath/$outDir/Data/Intermediate-files/DataTransformations/*.pdf $myPath/$outDir/Plots/
mv $myPath/$outDir/Plots/tsRNA.cond1-vs-cond2.high-distribution-score.pdf $myPath/$outDir/Plots/${condition1}_vs_${condition2}_tsRNAs.high-distribution-score.pdf
mv $myPath/$outDir/Plots/snomiRNA.cond1-vs-cond2.high-distribution-score.pdf $myPath/$outDir/Plots/${condition1}_vs_${condition2}_snomiRNAs.high-distribution-score.pdf

### Get mean and standard deviation
Rscript bin/Mean_Stdev.R $myPath/$outDir/Data/Intermediate-files/tsRNA.${condition1}_concatenated.depth $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.tsRNA.depth &
Rscript bin/Mean_Stdev.R $myPath/$outDir/Data/Intermediate-files/tsRNA.${condition2}_concatenated.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.tsRNA.depth &
Rscript bin/Mean_Stdev.R $myPath/$outDir/Data/Intermediate-files/snomiRNA.${condition1}_concatenated.depth $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth &
Rscript bin/Mean_Stdev.R $myPath/$outDir/Data/Intermediate-files/snomiRNA.${condition2}_concatenated.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth &
wait

### Get names of differentially expressed features
cat $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/*regulated.csv | grep -v ^,base | awk -F ',' '{print $1}' | awk -F ' ' '{print $1}' > $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt
grep ENS $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt > $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only_ENSGs.txt
grep -v ENS $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt | awk -F '-' '{print $2}' | sort | uniq > $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only_tRNAs.txt 
cat $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only_ENSGs.txt $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only_tRNAs.txt > $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only_short-names.txt
### Get names of features with highest distribution scores 
cat $myPath/$outDir/Data/Intermediate-files/DataTransformations/*high-distribution-score.sorted.txt | grep -v ^feat | awk '{print $1}' > $myPath/$outDir/Data/Intermediate-files/DataTransformations/High-distribution-scores_feature-names.txt 
### Get names of features that are likely cleaved
cat $myPath/$outDir/Data/Intermediate-files/*high-cleavage-score.txt | grep -v ^feat | awk '{print $1}' > $myPath/$outDir/Data/Intermediate-files/Potentially-cleaved-features_feature-names.txt

### Get unique set of differentially expressed features, features with high distribution scores, and high cleavage scores...
echo -e "#This is a collection of features that are differentially expressed, have large differences in distribution between the conditions, or are likely cleaved. Ordered alphanumerically." > $myPath/$outDir/Data/All-Features-Identified.txt
cat $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only_short-names.txt $myPath/$outDir/Data/Intermediate-files/DataTransformations/High-distribution-scores_feature-names.txt $myPath/$outDir/Data/Intermediate-files/Potentially-cleaved-features_feature-names.txt | sort | uniq >> $myPath/$outDir/Data/All-Features-Identified.txt

string_padder "Plotting any features identified in the analysis..."
### If there are differentially expressed/high distribution features, plot these:
if [[ $(wc -l < $myPath/$outDir/Data/All-Features-Identified.txt) -ge 2 ]]; then
	### Plot DEGs arg1 and 2 are inputs, arg 3 is list of differentially expressed genes, arg 4 is output pdf, 
	### arg 5 is mean coverage cutoff (plot features with coverage above this), arg 5 is GTF file for snomiRNAs (arg 5 is not given to tsRNA data)
	if [[ $Plots == "yes" ]]; then
		### Plot all tsRNAs and sno/miRNAs if -A 'yes' parameter selected
		Rscript bin/Bedgraph_plotter_DEGs.R $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.tsRNA.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.tsRNA.depth $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}_Features_All-tsRNAs.pdf 0 $Plots &
		Rscript bin/Bedgraph_plotter_DEGs.R $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}_Features_All-snomiRNAs.pdf 0 $Plots $snomiRNAGTF &
	fi
	# tsRNAs
	Rscript bin/Bedgraph_plotter_DEGs.R $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.tsRNA.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.tsRNA.depth $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}_Features_Differentially-expressed-tsRNAs.pdf 0 "no" &
	Rscript bin/Bedgraph_plotter_DEGs.R $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.tsRNA.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.tsRNA.depth $myPath/$outDir/Data/Intermediate-files/DataTransformations/High-distribution-scores_feature-names.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}_Features_High-distribution-tsRNAs.pdf 0 "no" &
	Rscript bin/Bedgraph_plotter_DEGs.R $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.tsRNA.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.tsRNA.depth $myPath/$outDir/Data/Intermediate-files/Potentially-cleaved-features_feature-names.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}_Features_Potentially-cleaved-tsRNAs.pdf 0 "no" &
	# snomiRNAs
	Rscript bin/Bedgraph_plotter_DEGs.R $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}_Features_Differentially-expressed-snomiRNAs.pdf 0 "no" $snomiRNAGTF &
	Rscript bin/Bedgraph_plotter_DEGs.R $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Data/Intermediate-files/DataTransformations/High-distribution-scores_feature-names.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}_Features_High-distribution-snomiRNAs.pdf 0 "no" $snomiRNAGTF &
	Rscript bin/Bedgraph_plotter_DEGs.R $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Data/Intermediate-files/Potentially-cleaved-features_feature-names.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}_Features_Potentially-cleaved-snomiRNAs.pdf 0 "no" $snomiRNAGTF &
	# Venn diagram
	Rscript bin/VennDiagram.R $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only_short-names.txt $myPath/$outDir/Data/Intermediate-files/DataTransformations/High-distribution-scores_feature-names.txt $myPath/$outDir/Data/Intermediate-files/Potentially-cleaved-features_feature-names.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}
	mv $myPath/$outDir/Plots/${condition1}_vs_${condition2}.intersect*txt $myPath/$outDir/Data/
else
	string_padder "No features of interest were identified."
	echo "There were no features identified. No plots were generated." >> $myPath/$outDir/Plots/${condition1}_vs_${condition2}_no-features-to-plot.txt
fi
wait

string_padder "Gathering RPM count files and cleaning up..."
### Gather RPM count files
awk '{print $1}' $outDir/Data/$filename.all-features.rpm.count > $outDir/Data/Intermediate-files/FCount.rpm.all-features # Get feature names
for f in $outDir/Data/*rpm.count; do
	awk '{print $2}' $f | paste $outDir/Data/Intermediate-files/FCount.rpm.all-features - >> $outDir/Data/Intermediate-files/FCount.rpm.temp
	mv $outDir/Data/Intermediate-files/FCount.rpm.temp $outDir/Data/Intermediate-files/FCount.rpm.all-features
done
rm $outDir/Data/*rpm.count
mv $outDir/Data/Intermediate-files/FCount.rpm.all-features $outDir/Data/FCount.all-features.rpm.count

### Move DESeq results to Data directory
mv $myPath/$outDir/Data/Intermediate-files/DE_Results/ $myPath/$outDir/Data/
cp $myPath/$outDir/Data/DE_Results/*pdf $myPath/$outDir/Plots/ #Copy PDFs to Plots dir
cp $myPath/$outDir/Data/DE_Results/DESeq2/*regulated.csv $myPath/$outDir/Data/ #Copy DESeq2 results to Data dir
mv $myPath/$outDir/Plots/*log $myPath/$outDir/Data/Intermediate-files/ #Move Venn Diagram log file to Intermediate-files dir 

### If -A parameter was provided, copy all plots to Plots dir
if [[ $Plots == "yes" ]]; then 
	mkdir $myPath/$outDir/Plots/Individual-Runs
	cp $myPath/$outDir/*/Data_and_Plots/*.pdf $myPath/$outDir/Plots/Individual-Runs/
fi

finished="Finished project analysis on $(date)"
string_padder "$finished"
