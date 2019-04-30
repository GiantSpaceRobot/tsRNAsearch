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
	-e	Optional (but recommended) CSV file containing file names and file groups (see examples in additional-files/)
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
fi

### Create dir substructure
for f in $inDir/*; do
	mkdir -p $outDir
	mkdir -p $outDir/Data
	mkdir -p $outDir/Data/Intermediate-files
	mkdir -p $outDir/Plots
	file_base=$(basename $f)
	filename="$( cut -d '.' -f 1 <<< "$file_base" )" 
	analysis="Beginning analysis of $filename using tsRNAsearch"
	string_padder $analysis
	tsRNAsearch.sh -g "$genome" -s "$f" -o "$outDir/$filename" -t "$threads" -A "$Plots" #>> "$outDir"/"$filename"_tiRNApipeline.log
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


### Determine if experiment layout file was provided or not. If not, try and figure out which files group together using R.
if [ ! "$expFile" ]; then
    string_padder "No experiment layout plan provided. This will now be created prior to the formal DESeq2 analysis."
    Rscript --vanilla bin/DESeq2_tsRNAsearch.R "$myPath/$outDir/Data/Intermediate-files/"
    condition1=$( awk -F ',' 'NR == 1 {print $2}' $myPath/$outDir/Data/Intermediate-files/predicted_exp_layout.csv ) # Get element in first row second column (condition)
    grep $condition1 $myPath/$outDir/Data/Intermediate-files/predicted_exp_layout.csv > $myPath/$outDir/Data/Intermediate-files/predicted_exp_layout_cond1.csv
    condition2=$( grep -v $condition1 $myPath/$outDir/Data/Intermediate-files/predicted_exp_layout.csv | awk -F ',' 'NR == 1 {print $2}') # Get the second condition using the inverse of the first one
    grep $condition2 $myPath/$outDir/Data/Intermediate-files/predicted_exp_layout.csv > $myPath/$outDir/Data/Intermediate-files/predicted_exp_layout_cond2.csv
	expFile="$myPath/$outDir/Data/Intermediate-files/predicted_exp_layout.csv"
else
    string_padder "An experiment layout plan was provided. Carrying out DESeq2 analysis now."
    Rscript --vanilla bin/DESeq2_tsRNAsearch.R "$expFile" "$myPath/$outDir/Data/Intermediate-files/"
    condition1=$( awk -F ',' 'NR == 1 {print $2}' "$expFile" ) # Get element in first row second column (condition)
    grep $condition1 "$expFile" > $myPath/$outDir/Data/Intermediate-files/predicted_exp_layout_cond1.csv
    condition2=$( grep -v $condition1 "$expFile" | awk -F ',' 'NR == 1 {print $2}') # Get the second condition using the inverse of the first one
    grep $condition2 "$expFile" > $myPath/$outDir/Data/Intermediate-files/predicted_exp_layout_cond2.csv
fi

### Concatenate the depth files for the replicates of condition 1
string_padder "Concatenating samtools depth files for condition 1..."
count=0
for fname in $(cat $myPath/$outDir/Data/Intermediate-files/predicted_exp_layout_cond1.csv); do
	count=$((count + 1))
	cond1="$(cut -d ',' -f1 <<< $fname)"
	cond1base="$(cut -d '.' -f1 <<< $cond1)"
	echo "Finding file for $fname ($cond1base)..."
	find $myPath/$outDir/ -type f -name "$cond1base*tiRNA.depth" -exec cp {} $myPath/$outDir/Data/Intermediate-files/condition1_file$count.tiRNA.depth \; & # Gather tiRNAs 
	find $myPath/$outDir/ -type f -name "$cond1base*snomiRNA.depth" -exec cp {} $myPath/$outDir/Data/Intermediate-files/condition1_file$count.snomiRNA.depth \; & # Gather sno/miRNAs
	wait
	mapped=$(grep "mapped" $myPath/$outDir/$cond1base/Stats.log | awk '{print $3}')
	echo "Converting raw counts to RPM..."
	python bin/Depth-to-Depth_RPM.py "$myPath/$outDir/Data/Intermediate-files/condition1_file$count.tiRNA.depth" "$mapped" "$myPath/$outDir/Data/Intermediate-files/condition1_file$count.tiRNA.depth.readspermil" &
	python bin/Depth-to-Depth_RPM.py "$myPath/$outDir/Data/Intermediate-files/condition1_file$count.snomiRNA.depth" "$mapped" "$myPath/$outDir/Data/Intermediate-files/condition1_file$count.snomiRNA.depth.readspermil" &
done
wait
echo "Transforming depth data to horizontal and vertical formats..."
### Concatenate data horizontally
paste $myPath/$outDir/Data/Intermediate-files/condition1_file*.tiRNA.depth.readspermil > $myPath/$outDir/Data/Intermediate-files/tiRNA.condition1_concatenated.depth &
paste $myPath/$outDir/Data/Intermediate-files/condition1_file*.snomiRNA.depth.readspermil > $myPath/$outDir/Data/Intermediate-files/snomiRNA.condition1_concatenated.depth &
### Concatenate data vertically
cat $myPath/$outDir/Data/Intermediate-files/condition1_file*.tiRNA.depth.readspermil | sort -k1,1 -k2,2n > $myPath/$outDir/Data/Intermediate-files/tiRNA.condition1_concatenated.depthVert &
cat $myPath/$outDir/Data/Intermediate-files/condition1_file*.snomiRNA.depth.readspermil | sort -k1,1 -k2,2n > $myPath/$outDir/Data/Intermediate-files/snomiRNA.condition1_concatenated.depthVert &
wait
 
### Concatenate the depth files for the replicates of condition 2 
string_padder "Concatenating samtools depth files for condition 2..."
count=0
for fname in $(cat $myPath/$outDir/Data/Intermediate-files/predicted_exp_layout_cond2.csv); do
	count=$((count + 1))
	cond2="$(cut -d ',' -f1 <<< $fname)"
	cond2base="$(cut -d '.' -f1 <<< $cond2)"
	echo "Finding file for $fname ($cond2base)..."
	find $myPath/$outDir/ -type f -name "$cond2base*tiRNA.depth" -exec cp {} $myPath/$outDir/Data/Intermediate-files/condition2_file$count.tiRNA.depth \; & # Gather tiRNAs
	find $myPath/$outDir/ -type f -name "$cond2base*snomiRNA.depth" -exec cp {} $myPath/$outDir/Data/Intermediate-files/condition2_file$count.snomiRNA.depth \; & # Gather snomiRNAs
	wait
	mapped=$(grep "mapped" $myPath/$outDir/$cond2base/Stats.log | awk '{print $3}')
	echo "Converting raw counts to RPM..."
	python bin/Depth-to-Depth_RPM.py "$myPath/$outDir/Data/Intermediate-files/condition2_file$count.tiRNA.depth" "$mapped" "$myPath/$outDir/Data/Intermediate-files/condition2_file$count.tiRNA.depth.readspermil" &
	python bin/Depth-to-Depth_RPM.py "$myPath/$outDir/Data/Intermediate-files/condition2_file$count.snomiRNA.depth" "$mapped" "$myPath/$outDir/Data/Intermediate-files/condition2_file$count.snomiRNA.depth.readspermil" &
done
wait
echo "Transforming depth data to horizontal and vertical formats..."
### Concatenate data horizontally
paste $myPath/$outDir/Data/Intermediate-files/condition2_file*.tiRNA.depth.readspermil > $myPath/$outDir/Data/Intermediate-files/tiRNA.condition2_concatenated.depth &
paste $myPath/$outDir/Data/Intermediate-files/condition2_file*.snomiRNA.depth.readspermil > $myPath/$outDir/Data/Intermediate-files/snomiRNA.condition2_concatenated.depth &
### Concatenate data vertically
cat $myPath/$outDir/Data/Intermediate-files/condition2_file*.tiRNA.depth.readspermil | sort -k1,1 -k2,2n > $myPath/$outDir/Data/Intermediate-files/tiRNA.condition2_concatenated.depthVert &
cat $myPath/$outDir/Data/Intermediate-files/condition2_file*.snomiRNA.depth.readspermil | sort -k1,1 -k2,2n > $myPath/$outDir/Data/Intermediate-files/snomiRNA.condition2_concatenated.depthVert &
wait

### Get distribution scores (standard deviation of RPM difference between samples multiplied 
### by standard deviation of percent difference between samples, divided by 1000) of the features
mkdir -p $myPath/$outDir/Data/Intermediate-files/DataTransformations

string_padder "tiRNAs: Generating Distribution Scores..."
### tiRNAs
### Calculate mean
python bin/MeanCalculator.py $myPath/$outDir/Data/Intermediate-files/tiRNA.condition1_concatenated.depth $myPath/$outDir/Data/Intermediate-files/DataTransformations/tiRNA.condition1_concatenated.depth.mean & 
python bin/MeanCalculator.py $myPath/$outDir/Data/Intermediate-files/tiRNA.condition2_concatenated.depth $myPath/$outDir/Data/Intermediate-files/DataTransformations/tiRNA.condition2_concatenated.depth.mean &
wait
### Sort output
sort -k1,1 -k2,2n $myPath/$outDir/Data/Intermediate-files/DataTransformations/tiRNA.condition1_concatenated.depth.mean > $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_tiRNA.condition1_concatenated.depth.mean &
sort -k1,1 -k2,2n $myPath/$outDir/Data/Intermediate-files/DataTransformations/tiRNA.condition2_concatenated.depth.mean > $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_tiRNA.condition2_concatenated.depth.mean &
wait
### Concatenate horizontally to make a dataframe
paste $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_tiRNA.condition1_concatenated.depth.mean $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_tiRNA.condition2_concatenated.depth.mean > $myPath/$outDir/Data/Intermediate-files/DataTransformations/tiRNA.cond1-vs-cond2.mean
### Calculate StdDev
python bin/Mean-to-RelativeDifference.py $myPath/$outDir/Data/Intermediate-files/DataTransformations/tiRNA.cond1-vs-cond2.mean $myPath/$outDir/Data/Intermediate-files/DataTransformations/tiRNA.cond1-vs-cond2.stddev
### Calculate distribution scores
Rscript bin/DistributionScore.R $myPath/$outDir/Data/Intermediate-files/DataTransformations/tiRNA.cond1-vs-cond2.stddev $myPath/$outDir/Data/Intermediate-files/DataTransformations/tiRNA.cond1-vs-cond2 $snomiRNAGTF
### Sort the output but not the header
cat $myPath/$outDir/Data/Intermediate-files/DataTransformations/tiRNA.cond1-vs-cond2.high-distribution-score.txt | awk 'NR<2{print $0;next}{print $0| "sort -k5,5nr"}' > $myPath/$outDir/Data/Intermediate-files/DataTransformations/tiRNA.cond1-vs-cond2.high-distribution-score.sorted.txt
cat $myPath/$outDir/Data/Intermediate-files/DataTransformations/tiRNA.cond1-vs-cond2.all-features.txt | awk 'NR<2{print $0;next}{print $0| "sort -k11,11nr"}' > $myPath/$outDir/Data/Intermediate-files/DataTransformations/tiRNA.cond1-vs-cond2.all-features.sorted.txt

string_padder "snoRNAs/miRNAs: Generating Distribution Scores..."
### snomiRNAs
### Calculate mean
python bin/MeanCalculator.py $myPath/$outDir/Data/Intermediate-files/snomiRNA.condition1_concatenated.depth $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.condition1_concatenated.depth.mean &
python bin/MeanCalculator.py $myPath/$outDir/Data/Intermediate-files/snomiRNA.condition2_concatenated.depth $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.condition2_concatenated.depth.mean &
wait
### Sort output
sort -k1,1 -k2,2n $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.condition1_concatenated.depth.mean > $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_snomiRNA.condition1_concatenated.depth.mean &
sort -k1,1 -k2,2n $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.condition2_concatenated.depth.mean > $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_snomiRNA.condition2_concatenated.depth.mean &
wait
### Concatenate horizontally to make a dataframe
paste $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_snomiRNA.condition1_concatenated.depth.mean $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_snomiRNA.condition2_concatenated.depth.mean > $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.mean
### Calculate StdDev
python bin/Mean-to-RelativeDifference.py $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.mean $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.stddev
### Calculate distribution score
Rscript bin/DistributionScore.R $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.stddev $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2 $snomiRNAGTF
### Sort the output but not the header
cat $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.high-distribution-score.txt | awk 'NR<2{print $0;next}{print $0| "sort -k6,6nr"}' > $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.high-distribution-score.sorted.txt
#cat $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.all-features.txt | awk 'NR<2{print $0;next}{print $0| "sort -k14,14nr"}' > $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.all-features.sorted.txt

### Copy txt files illustrating features with very different distributions to the Results directory:
cp $myPath/$outDir/Data/Intermediate-files/DataTransformations/tiRNA.cond1-vs-cond2.high-distribution-score.sorted.txt $myPath/$outDir/Data/${condition1}_vs_${condition2}_High-distribution-tiRNAs.txt 
cp $myPath/$outDir/Data/Intermediate-files/DataTransformations/snomiRNA.cond1-vs-cond2.high-distribution-score.sorted.txt $myPath/$outDir/Data/${condition1}_vs_${condition2}_High-distribution-snomiRNAs.txt


string_padder "Generating Cleavage Scores..." 
### Test whether certain features are cleaved in one condition vs the other
Rscript bin/CleavageScore.R $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_tiRNA.condition1_concatenated.depth.mean $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_tiRNA.condition2_concatenated.depth.mean $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_tiRNAs $snomiRNAGTF &
Rscript bin/CleavageScore.R $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_snomiRNA.condition1_concatenated.depth.mean $myPath/$outDir/Data/Intermediate-files/DataTransformations/sorted_snomiRNA.condition2_concatenated.depth.mean $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_snomiRNAs $snomiRNAGTF &
wait
### Sort the output but not the header
cat $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_tiRNAs.high-cleavage-score.txt | awk 'NR<2{print $0;next}{print $0| "sort -k9,9nr"}' > $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_High-cleavage-tiRNAs.txt
cat $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_snomiRNAs.high-cleavage-score.txt | awk 'NR<2{print $0;next}{print $0| "sort -k10,10nr"}' > $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_High-cleavage-snomiRNAs.txt
### Copy sorted high distribution files to the Data dir
cp $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_High-cleavage-snomiRNAs.txt $myPath/$outDir/Data/
cp $myPath/$outDir/Data/Intermediate-files/${condition1}_vs_${condition2}_High-cleavage-tiRNAs.txt $myPath/$outDir/Data/
### Copy PDF files to Plots dir
cp $myPath/$outDir/Data/Intermediate-files/*.pdf $myPath/$outDir/Plots/
cp $myPath/$outDir/Data/Intermediate-files/DataTransformations/*.pdf $myPath/$outDir/Plots/
mv $myPath/$outDir/Plots/tiRNA.cond1-vs-cond2.high-distribution-score.pdf $myPath/$outDir/Plots/${condition1}_vs_${condition2}_tiRNAs.high-distribution-score.pdf
mv $myPath/$outDir/Plots/snomiRNA.cond1-vs-cond2.high-distribution-score.pdf $myPath/$outDir/Plots/${condition1}_vs_${condition2}_snomiRNAs.high-distribution-score.pdf

### Get mean and standard deviation
Rscript bin/Mean_Stdev.R $myPath/$outDir/Data/Intermediate-files/tiRNA.condition1_concatenated.depth $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.tiRNA.depth &
Rscript bin/Mean_Stdev.R $myPath/$outDir/Data/Intermediate-files/tiRNA.condition2_concatenated.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.tiRNA.depth &
Rscript bin/Mean_Stdev.R $myPath/$outDir/Data/Intermediate-files/snomiRNA.condition1_concatenated.depth $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth &
Rscript bin/Mean_Stdev.R $myPath/$outDir/Data/Intermediate-files/snomiRNA.condition2_concatenated.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth &
wait

### Get names of differentially expressed features
cat $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/*regulated.csv | grep -v ^,base | awk -F ',' '{print $1}' > $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt
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
	### arg 5 is mean coverage cutoff (plot features with coverage above this), arg 5 is GTF file for snomiRNAs (arg 5 is not given to tiRNA data)
	# tiRNAs
	Rscript bin/Bedgraph_plotter_DEGs.R $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.tiRNA.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.tiRNA.depth $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}_Features_Differentially-expressed-tiRNAs.pdf 0 &
	Rscript bin/Bedgraph_plotter_DEGs.R $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.tiRNA.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.tiRNA.depth $myPath/$outDir/Data/Intermediate-files/DataTransformations/High-distribution-scores_feature-names.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}_Features_High-distribution-tiRNAs.pdf 0 &
	Rscript bin/Bedgraph_plotter_DEGs.R $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.tiRNA.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.tiRNA.depth $myPath/$outDir/Data/Intermediate-files/Potentially-cleaved-features_feature-names.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}_Features_Potentially-cleaved-tiRNAs.pdf 0 &
	# snomiRNAs
	Rscript bin/Bedgraph_plotter_DEGs.R $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}_Features_Differentially-expressed-snomiRNAs.pdf 0 $snomiRNAGTF &
	Rscript bin/Bedgraph_plotter_DEGs.R $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Data/Intermediate-files/DataTransformations/High-distribution-scores_feature-names.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}_Features_High-distribution-snomiRNAs.pdf 0 $snomiRNAGTF &
	Rscript bin/Bedgraph_plotter_DEGs.R $myPath/$outDir/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Data/Intermediate-files/Potentially-cleaved-features_feature-names.txt $myPath/$outDir/Plots/${condition1}_vs_${condition2}_Features_Potentially-cleaved-snomiRNAs.pdf 0 $snomiRNAGTF &
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
cp $myPath/$outDir/Data/DE_Results/*pdf $myPath/$outDir/Plots/

finished="Finished project analysis on $(date)"
string_padder "$finished"
