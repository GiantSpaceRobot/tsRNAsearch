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
	-e	Optional (but recommended) CSV file containing file names and file groups (see examples in ./additional-files/)
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
if [ "$expFile" ]; then
    if [[ ! $expFile = /* ]]; then
        echo "Error: File path must absolute. Please specify the full path for the experiment layout file."
        exit 1
    fi
fi


##### Start of pipeline #####
echo "Started project analysis on $(date)"
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

### If -A parameter was not provided, default is to only plot differentially expressed features
if [ ! "$Plots" ]; then
	Plots="no"
fi

### Create dir substructure
for f in $inDir/*; do
	mkdir -p $outDir
	mkdir -p $outDir/Results
	mkdir -p $outDir/Results/Data
	mkdir -p $outDir/Results/Data/Intermediate-files
	mkdir -p $outDir/Results/Plots
	file_base=$(basename $f)
	filename="$( cut -d '.' -f 1 <<< "$file_base" )" 
	#if [ "$tophat" ]; then # Use tophat2
	#	./tiRNA-pipeline.sh -s "$f" -o "$outDir/Results/$filename" -p "$threads" -T -A "$Plots" #>> "$outDir"/"$filename"_tiRNApipeline.log
	#else # Use HISAT2
	./tsRNAsearch.sh -g "$genome" -s "$f" -o "$outDir/Results/$filename" -p "$threads" -A "$Plots" #>> "$outDir"/"$filename"_tiRNApipeline.log
	#fi
	wait
	#cp $outDir/Results/$filename/Data_and_Plots/*pdf $outDir/Results/Plots/
	cat $outDir/Results/$filename/FCount-count-output/*.count | grep -v ^__ | sort -k1,1 > $outDir/Results/Data/Intermediate-files/$filename.all-features.count	
	readsMapped=$(awk '{sum+=$2} END{print sum;}' $outDir/Results/Data/Intermediate-files/$filename.all-features.count)
	cp $outDir/Results/$filename/Data_and_Plots/$filename.all-features.rpm.count $outDir/Results/Data/	
	cp $outDir/Results/$filename/FCount-count-output/$filename.all-features.count $outDir/Results/Data/Intermediate-files/ 
done

### Gather raw count files
awk '{print $1}' $outDir/Results/Data/Intermediate-files/$filename.all-features.count > $outDir/Results/Data/Intermediate-files/FCount.all-features # Get feature names
for f in $outDir/Results/Data/Intermediate-files/*count; do
	awk '{print $2}' $f | paste $outDir/Results/Data/Intermediate-files/FCount.all-features - >> $outDir/Results/Data/Intermediate-files/FCount.temp
	mv $outDir/Results/Data/Intermediate-files/FCount.temp $outDir/Results/Data/Intermediate-files/FCount.all-features
done
mv $outDir/Results/Data/Intermediate-files/FCount.all-features $outDir/Results/Data/FCount.all-features.raw.count


### Determine if experiment layout file was provided or not. If not, try and figure out which files group together using R.
myPath=$(pwd) #Get working dir to recreate full path for R script execution
if [ ! "$expFile" ]; then
    echo "No experiment layout plan provided. This will now be created prior to the formal DESeq2 analysis."
    Rscript --vanilla scripts/DESeq2_tiRNA-pipeline-v2.R "$myPath/$outDir/Results/Data/Intermediate-files/"
    condition1=$( awk -F ',' 'NR == 1 {print $2}' $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout.csv ) # Get element in first row second column (condition)
    grep $condition1 $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout.csv > $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout_cond1.csv
    condition2=$( grep -v $condition1 $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout.csv | awk -F ',' 'NR == 1 {print $2}') # Get the second condition using the inverse of the first one
    grep $condition2 $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout.csv > $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout_cond2.csv
	expFile="$myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout.csv"
else
    echo "An experiment layout plan was provided. Carrying out DESeq2 analysis now."
    Rscript --vanilla scripts/DESeq2_tiRNA-pipeline-v2.R "$expFile" "$myPath/$outDir/Results/Data/Intermediate-files/"
    condition1=$( awk -F ',' 'NR == 1 {print $2}' "$expFile" ) # Get element in first row second column (condition)
    grep $condition1 "$expFile" > $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout_cond1.csv
    condition2=$( grep -v $condition1 "$expFile" | awk -F ',' 'NR == 1 {print $2}') # Get the second condition using the inverse of the first one
    grep $condition2 "$expFile" > $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout_cond2.csv
fi

### Concatenate the depth files for the replicates of condition 1
count=0
for fname in $(cat $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout_cond1.csv); do
	count=$((count + 1))
	cond1="$(cut -d ',' -f1 <<< $fname)"
	cond1base="$(cut -d '.' -f1 <<< $cond1)"
	find $myPath/$outDir/ -type f -name "$cond1base*tiRNA.depth" -exec cp {} $myPath/$outDir/Results/Data/Intermediate-files/condition1_file$count.tiRNA.depth \; & # Gather tiRNAs 
	find $myPath/$outDir/ -type f -name "$cond1base*snomiRNA.depth" -exec cp {} $myPath/$outDir/Results/Data/Intermediate-files/condition1_file$count.snomiRNA.depth \; & # Gather sno/miRNAs
	wait
	mapped=$(grep "mapped" $myPath/$outDir/Results/$cond1base/Stats.log | awk '{print $3}')
	python scripts/Depth-to-Depth_RPM.py "$myPath/$outDir/Results/Data/Intermediate-files/condition1_file$count.tiRNA.depth" "$mapped" "$myPath/$outDir/Results/Data/Intermediate-files/condition1_file$count.tiRNA.depth.readspermil" &
	python scripts/Depth-to-Depth_RPM.py "$myPath/$outDir/Results/Data/Intermediate-files/condition1_file$count.snomiRNA.depth" "$mapped" "$myPath/$outDir/Results/Data/Intermediate-files/condition1_file$count.snomiRNA.depth.readspermil" &
done
wait
### Concatenate data horizontally
paste $myPath/$outDir/Results/Data/Intermediate-files/condition1_file*.tiRNA.depth.readspermil > $myPath/$outDir/Results/Data/Intermediate-files/tiRNA.condition1_concatenated.depth &
paste $myPath/$outDir/Results/Data/Intermediate-files/condition1_file*.snomiRNA.depth.readspermil > $myPath/$outDir/Results/Data/Intermediate-files/snomiRNA.condition1_concatenated.depth &
### Concatenate data vertically
cat $myPath/$outDir/Results/Data/Intermediate-files/condition1_file*.tiRNA.depth.readspermil | sort -k1,1 -k2,2n > $myPath/$outDir/Results/Data/Intermediate-files/tiRNA.condition1_concatenated.depthVert &
cat $myPath/$outDir/Results/Data/Intermediate-files/condition1_file*.snomiRNA.depth.readspermil | sort -k1,1 -k2,2n > $myPath/$outDir/Results/Data/Intermediate-files/snomiRNA.condition1_concatenated.depthVert &
wait
 
### Concatenate the depth files for the replicates of condition 2 
count=0
for fname in $(cat $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout_cond2.csv); do
	count=$((count + 1))
	cond2="$(cut -d ',' -f1 <<< $fname)"
	cond2base="$(cut -d '.' -f1 <<< $cond2)"
	find $myPath/$outDir/ -type f -name "$cond2base*tiRNA.depth" -exec cp {} $myPath/$outDir/Results/Data/Intermediate-files/condition2_file$count.tiRNA.depth \; & # Gather tiRNAs
	find $myPath/$outDir/ -type f -name "$cond2base*snomiRNA.depth" -exec cp {} $myPath/$outDir/Results/Data/Intermediate-files/condition2_file$count.snomiRNA.depth \; & # Gather snomiRNAs
	wait
	mapped=$(grep "mapped" $myPath/$outDir/Results/$cond2base/Stats.log | awk '{print $3}')
	python scripts/Depth-to-Depth_RPM.py "$myPath/$outDir/Results/Data/Intermediate-files/condition2_file$count.tiRNA.depth" "$mapped" "$myPath/$outDir/Results/Data/Intermediate-files/condition2_file$count.tiRNA.depth.readspermil" &
	python scripts/Depth-to-Depth_RPM.py "$myPath/$outDir/Results/Data/Intermediate-files/condition2_file$count.snomiRNA.depth" "$mapped" "$myPath/$outDir/Results/Data/Intermediate-files/condition2_file$count.snomiRNA.depth.readspermil" &
done
wait
### Concatenate data horizontally
paste $myPath/$outDir/Results/Data/Intermediate-files/condition2_file*.tiRNA.depth.readspermil > $myPath/$outDir/Results/Data/Intermediate-files/tiRNA.condition2_concatenated.depth &
paste $myPath/$outDir/Results/Data/Intermediate-files/condition2_file*.snomiRNA.depth.readspermil > $myPath/$outDir/Results/Data/Intermediate-files/snomiRNA.condition2_concatenated.depth &
### Concatenate data vertically
cat $myPath/$outDir/Results/Data/Intermediate-files/condition2_file*.tiRNA.depth.readspermil | sort -k1,1 -k2,2n > $myPath/$outDir/Results/Data/Intermediate-files/tiRNA.condition2_concatenated.depthVert &
cat $myPath/$outDir/Results/Data/Intermediate-files/condition2_file*.snomiRNA.depth.readspermil | sort -k1,1 -k2,2n > $myPath/$outDir/Results/Data/Intermediate-files/snomiRNA.condition2_concatenated.depthVert &
wait

### Get distribution scores (standard deviation of RPM difference between samples multiplied 
### by standard deviation of percent difference between samples, divided by 1000) of the features
mkdir -p $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score

### tiRNAs
### Calculate mean
python scripts/MeanCalculator.py $myPath/$outDir/Results/Data/Intermediate-files/tiRNA.condition1_concatenated.depth $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.condition1_concatenated.depth.mean & 
python scripts/MeanCalculator.py $myPath/$outDir/Results/Data/Intermediate-files/tiRNA.condition2_concatenated.depth $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.condition2_concatenated.depth.mean &
wait
### Sort output
sort -k1,1 -k2,2n $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.condition1_concatenated.depth.mean > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/sorted_tiRNA.condition1_concatenated.depth.mean &
sort -k1,1 -k2,2n $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.condition2_concatenated.depth.mean > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/sorted_tiRNA.condition2_concatenated.depth.mean &
wait
### Concatenate horizontally to make a dataframe
paste $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/sorted_tiRNA.condition1_concatenated.depth.mean $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/sorted_tiRNA.condition2_concatenated.depth.mean > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.mean
### Calculate StdDev
python scripts/Mean-to-RelativeDifference.py $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.mean $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.stddev
### Calculate distribution scores
Rscript scripts/Distribution-score-v2.R $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.stddev $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2
### Sort the output but not the header
cat $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.different-distributions.txt | awk 'NR<2{print $0;next}{print $0| "sort -k11,11nr"}' > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.different-distributions.sorted.txt
cat $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.all-features.txt | awk 'NR<2{print $0;next}{print $0| "sort -k11,11nr"}' > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.all-features.sorted.txt

### snomiRNAs
### Calculate mean
python scripts/MeanCalculator.py $myPath/$outDir/Results/Data/Intermediate-files/snomiRNA.condition1_concatenated.depth $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.condition1_concatenated.depth.mean &
python scripts/MeanCalculator.py $myPath/$outDir/Results/Data/Intermediate-files/snomiRNA.condition2_concatenated.depth $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.condition2_concatenated.depth.mean &
wait
### Sort output
sort -k1,1 -k2,2n $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.condition1_concatenated.depth.mean > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/sorted_snomiRNA.condition1_concatenated.depth.mean &
sort -k1,1 -k2,2n $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.condition2_concatenated.depth.mean > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/sorted_snomiRNA.condition2_concatenated.depth.mean &
wait
### Concatenate horizontally to make a dataframe
paste $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/sorted_snomiRNA.condition1_concatenated.depth.mean $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/sorted_snomiRNA.condition2_concatenated.depth.mean > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2.mean
### Calculate StdDev
python scripts/Mean-to-RelativeDifference.py $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2.mean $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2.stddev
### Calculate distribution score
Rscript scripts/Distribution-score-v2.R $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2.stddev $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2
### Sort the output but not the header
cat $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2.different-distributions.txt | awk 'NR<2{print $0;next}{print $0| "sort -k11,11nr"}' > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2.different-distributions.sorted.txt
cat $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2.all-features.txt | awk 'NR<2{print $0;next}{print $0| "sort -k11,11nr"}' > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2.all-features.sorted.txt

### Copy txt files illustrating features with very different distributions to the Results directory:
cp $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.different-distributions.sorted.txt $myPath/$outDir/Results/Data/${condition1}_vs_${condition2}_High-distribution-tiRNAs.txt 
cp $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2.different-distributions.sorted.txt $myPath/$outDir/Results/Data/${condition1}_vs_${condition2}_High-distribution-snomiRNAs.txt

### Test whether certain features are cleaved in one condition vs the other
Rscript scripts/FeatureCleavageScore.R $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/sorted_tiRNA.condition1_concatenated.depth.mean $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/sorted_tiRNA.condition2_concatenated.depth.mean $myPath/$outDir/Results/Data/Intermediate-files/${condition1}_vs_${condition2}_tiRNAs &
Rscript scripts/FeatureCleavageScore.R $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/sorted_snomiRNA.condition1_concatenated.depth.mean $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/sorted_snomiRNA.condition2_concatenated.depth.mean $myPath/$outDir/Results/Data/Intermediate-files/${condition1}_vs_${condition2}_snomiRNAs &
wait
cp $myPath/$outDir/Results/Data/Intermediate-files/${condition1}_vs_${condition2}_tiRNAs.potentially-cleaved-features.txt $myPath/$outDir/Results/Data/
cp $myPath/$outDir/Results/Data/Intermediate-files/${condition1}_vs_${condition2}_snomiRNAs.potentially-cleaved-features.txt $myPath/$outDir/Results/Data/
cp $myPath/$outDir/Results/Data/Intermediate-files/${condition1}_vs_${condition2}_tiRNAs.potentially-cleaved-features.pdf $myPath/$outDir/Results/Plots/
cp $myPath/$outDir/Results/Data/Intermediate-files/${condition1}_vs_${condition2}_snomiRNAs.potentially-cleaved-features.pdf $myPath/$outDir/Results/Plots/

### Get mean and standard deviation
Rscript scripts/Mean_Stdev.R $myPath/$outDir/Results/Data/Intermediate-files/tiRNA.condition1_concatenated.depth $myPath/$outDir/Results/Data/Intermediate-files/condition1_concatenated_mean_stdev.tiRNA.depth &
Rscript scripts/Mean_Stdev.R $myPath/$outDir/Results/Data/Intermediate-files/tiRNA.condition2_concatenated.depth $myPath/$outDir/Results/Data/Intermediate-files/condition2_concatenated_mean_stdev.tiRNA.depth &
Rscript scripts/Mean_Stdev.R $myPath/$outDir/Results/Data/Intermediate-files/snomiRNA.condition1_concatenated.depth $myPath/$outDir/Results/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth &
Rscript scripts/Mean_Stdev.R $myPath/$outDir/Results/Data/Intermediate-files/snomiRNA.condition2_concatenated.depth $myPath/$outDir/Results/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth &
wait

### Get names of differentially expressed features
cat $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/*regulated.csv | grep -v ^,base | awk -F ',' '{print $1}' > $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt
grep ENS $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt > $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only_ENSGs.txt
grep -v ENS $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt | awk -F '-' '{print $2}' | sort | uniq > $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only_tRNAs.txt 
cat $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only_ENSGs.txt $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only_tRNAs.txt > $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only_short-names.txt
### Get names of features with highest distribution scores 
cat $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/*different-distributions.sorted.txt | grep -v ^feat | awk '{print $1}' > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/High-distribution-scores_feature-names.txt 
### Get names of features that are likely cleaved
cat $myPath/$outDir/Results/Data/Intermediate-files/*potentially-cleaved-features.txt | grep -v ^feat | awk '{print $1}' > $myPath/$outDir/Results/Data/Intermediate-files/Potentially-cleaved-features_feature-names.txt

### Get unique set of differentially expressed features and features with high distribution scores
echo -e "#This is a collection of features that are differentially expressed, have large differences in distribution between the conditions, or are likely cleaved. Ordered alphanumerically." > $myPath/$outDir/Results/Data/All-Features-Identified.txt
cat $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only_short-names.txt $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/High-distribution-scores_feature-names.txt $myPath/$outDir/Results/Data/Intermediate-files/Potentially-cleaved-features_feature-names.txt | sort | uniq >> $myPath/$outDir/Results/Data/All-Features-Identified.txt

### If there are differentially expressed/high distribution features, plot these:
if [[ $(wc -l < $myPath/$outDir/Results/Data/All-Features-Identified.txt) -ge 2 ]]; then
	### Plot DEGs arg1 and 2 are inputs, arg 3 is list of differentially expressed genes, arg 4 is output pdf, 
	### arg 5 is mean coverage cutoff (plot features with coverage above this), arg 5 is GTF file for snomiRNAs (arg 5 is not given to tiRNA data)
	# tiRNAs
	Rscript scripts/Bedgraph_plotter_DEGs-v4.R $myPath/$outDir/Results/Data/Intermediate-files/condition1_concatenated_mean_stdev.tiRNA.depth $myPath/$outDir/Results/Data/Intermediate-files/condition2_concatenated_mean_stdev.tiRNA.depth $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt $myPath/$outDir/Results/Plots/${condition1}_vs_${condition2}_Features_Differentially-expressed-tiRNAs.pdf 0 &
	Rscript scripts/Bedgraph_plotter_DEGs-v4.R $myPath/$outDir/Results/Data/Intermediate-files/condition1_concatenated_mean_stdev.tiRNA.depth $myPath/$outDir/Results/Data/Intermediate-files/condition2_concatenated_mean_stdev.tiRNA.depth $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/High-distribution-scores_feature-names.txt $myPath/$outDir/Results/Plots/${condition1}_vs_${condition2}_Features_High-distribution-tiRNAs.pdf 0 &
	Rscript scripts/Bedgraph_plotter_DEGs-v4.R $myPath/$outDir/Results/Data/Intermediate-files/condition1_concatenated_mean_stdev.tiRNA.depth $myPath/$outDir/Results/Data/Intermediate-files/condition2_concatenated_mean_stdev.tiRNA.depth $myPath/$outDir/Results/Data/Intermediate-files/Potentially-cleaved-features_feature-names.txt $myPath/$outDir/Results/Plots/${condition1}_vs_${condition2}_Features_Potentially-cleaved-tiRNAs.pdf 0 &
	# snomiRNAs
	Rscript scripts/Bedgraph_plotter_DEGs-v4.R $myPath/$outDir/Results/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Results/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt $myPath/$outDir/Results/Plots/${condition1}_vs_${condition2}_Features_Differentially-expressed-snomiRNAs.pdf 0 $snomiRNAGTF &
	Rscript scripts/Bedgraph_plotter_DEGs-v4.R $myPath/$outDir/Results/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Results/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/High-distribution-scores_feature-names.txt $myPath/$outDir/Results/Plots/${condition1}_vs_${condition2}_Features_High-distribution-snomiRNAs.pdf 0 $snomiRNAGTF &
	Rscript scripts/Bedgraph_plotter_DEGs-v4.R $myPath/$outDir/Results/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Results/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Results/Data/Intermediate-files/Potentially-cleaved-features_feature-names.txt $myPath/$outDir/Results/Plots/${condition1}_vs_${condition2}_Features_Potentially-cleaved-snomiRNAs.pdf 0 $snomiRNAGTF &
	# Venn diagram
	Rscript scripts/VennDiagram.R $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only_short-names.txt $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/High-distribution-scores_feature-names.txt $myPath/$outDir/Results/Data/Intermediate-files/Potentially-cleaved-features_feature-names.txt $myPath/$outDir/Results/Plots/${condition1}_vs_${condition2}
	mv $myPath/$outDir/Results/Plots/${condition1}_vs_${condition2}.intersect*txt $myPath/$outDir/Results/Data/
else
	echo "There were no differentially expressed features. No plots were generated." >> $myPath/$outDir/Results/Plots/${condition1}_vs_${condition2}_no-features-to-plot.txt
fi
wait

### Gather RPM count files
awk '{print $1}' $outDir/Results/Data/$filename.all-features.rpm.count > $outDir/Results/Data/Intermediate-files/FCount.rpm.all-features # Get feature names
for f in $outDir/Results/Data/*rpm.count; do
	awk '{print $2}' $f | paste $outDir/Results/Data/Intermediate-files/FCount.rpm.all-features - >> $outDir/Results/Data/Intermediate-files/FCount.rpm.temp
	mv $outDir/Results/Data/Intermediate-files/FCount.rpm.temp $outDir/Results/Data/Intermediate-files/FCount.rpm.all-features
done
rm $outDir/Results/Data/*rpm.count
mv $outDir/Results/Data/Intermediate-files/FCount.rpm.all-features $outDir/Results/Data/FCount.all-features.rpm.count

### Move DESeq results to Data directory
mv $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/ $myPath/$outDir/Results/Data/
cp $myPath/$outDir/Results/Data/DE_Results/*pdf $myPath/$outDir/Results/Plots/

echo "Finished project analysis on $(date)"

