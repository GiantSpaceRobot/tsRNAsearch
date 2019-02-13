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
	-A	Plot all features? yes/no {default: no (only plot differentially expressed features)}
	-T	Use Tophat2 instead of HISAT2 for alignment steps {default is HISAT2}

	" 1>&2; }

while getopts ":hTpA:t:d:e:o:" o; do
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
			;;
		e)
			expFile="$OPTARG"
			;;
		o)
			outDir="$OPTARG"
			;;
		t)
			CPUs="$OPTARG"
			;;
		A)
			Plots="$OPTARG"
			;;
		T)
			tophat="Yes"
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
#if [[ ! $outDir = /* ]]; then
#	echo "Error: File paths must absolute. Please specify the full path for the output directory."
#	exit 1
#fi

### If the pathname specified by $expFile does not begin with a slash, quit (we need full path name)
if [ "$expFile" ]; then
    if [[ ! $expFile = /* ]]; then
        echo "Error: File paths must absolute. Please specify the full path for the experiment layout file."
        exit 1
    fi
fi


##### Start of pipeline #####
echo "Started at $(date)"
StartTime="Pipeline initiated at $(date)"

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
	if [ "$tophat" ]; then # Use tophat2
		./tiRNA-pipeline.sh -s "$f" -o "$outDir/Results/$filename" -p "$CPUs" -T -A "$Plots" #>> "$outDir"/"$filename"_tiRNApipeline.log
	else # Use HISAT2
		./tiRNA-pipeline.sh -s "$f" -o "$outDir/Results/$filename" -p "$CPUs" -A "$Plots" #>> "$outDir"/"$filename"_tiRNApipeline.log
	fi
	wait
	cp $outDir/Results/$filename/Data_and_Plots/* $outDir/Results/Plots/
	#cat $outDir/Results/$filename/HTSeq-count-output/*.count | grep -v ^__ | sort -k1,1 > $outDir/Results/Data/Intermediate-files/$filename.all_features.count	
	#sed -i '1s/^/Features\t'"$filename"'\n/' $outDir/Results/Data/Intermediate-files/$filename.all_features.count # Add column headers
	readsMapped=$(awk '{sum+=$2} END{print sum;}' $outDir/Results/Data/Intermediate-files/$filename.all_features.count)
	#echo "Reads mapped: $readsMapped" >> $outDir/Results/$filename/Stats.log
	cp $outDir/Results/$filename/Data_and_Plots/$filename.all-features.rpm.count $outDir/Results/Data/Intermediate-files/$filename.all_features.count	
	
	
	#python scripts/Depth-to-RPM.py $outDir/Results/$filename/tRNA-alignment/accepted_hits_sorted.depth $readsMapped $outDir/Results/Data/Intermediate-files/$filename.depth
	
	
	### Normalise HTSeq-count files by total reads mapped to get RPM (reads per million) 
	#python scripts/HTSeq-to-RPM.py $outDir/Results/$filename/HTSeq-count-output/tRNA-alignment.count $readsMapped $outDir/Results/$filename/HTSeq-to-RPM/tRNA-alignment.RPM &
	#python scripts/HTSeq-to-RPM.py $outDir/Results/$filename/HTSeq-count-output/snomiRNA-alignment.count $readsMapped $outDir/Results/$filename/HTSeq-to-RPM/snomiRNA-alignment.RPM &
	#python scripts/HTSeq-to-RPM.py $outDir/Results/$filename/HTSeq-count-output/mRNA-ncRNA-alignment.count $readsMapped $outDir/Results/$filename/HTSeq-to-RPM/mRNA-ncRNA-alignment.RPM &
	#wait
done

### Gather count files
awk '{print $1}' $outDir/Results/Data/Intermediate-files/$filename.all_features.count > $outDir/Results/Data/Intermediate-files/HTSeq.all_features # Get feature names
for f in $outDir/Results/Data/Intermediate-files/*count; do
	awk '{print $2}' $f | paste $outDir/Results/Data/Intermediate-files/HTSeq.all_features - >> $outDir/Results/Data/Intermediate-files/HTSeq.temp
	mv $outDir/Results/Data/Intermediate-files/HTSeq.temp $outDir/Results/Data/Intermediate-files/HTSeq.all_features
done

mv $outDir/Results/Data/Intermediate-files/HTSeq.all_features $outDir/Results/Plots/HTSeq.all_features.count


### Determine if experiment layout file was provided or not. If not, try and figure out which files group together using R.
myPath=$(pwd) #Get working dir to recreate full path for R script execution
if [ ! "$expFile" ]; then
    echo "No experiment layout plan provided. This will now be created prior to the formal DESeq2 analysis."
    Rscript --vanilla scripts/DESeq2_tiRNA-pipeline.R "$myPath/$outDir/Results/Data/Intermediate-files/"
    condition1=$( awk -F ',' 'NR == 1 {print $2}' $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout.csv ) # Get element in first row second column (condition)
    grep $condition1 $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout.csv > $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout_cond1.csv
    condition2=$( grep -v $condition1 $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout.csv | awk -F ',' 'NR == 1 {print $2}') # Get the second condition using the inverse of the first one
    grep $condition2 $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout.csv > $myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout_cond2.csv
	expFile="$myPath/$outDir/Results/Data/Intermediate-files/predicted_exp_layout.csv"
else
    echo "An experiment layout plan was provided. Carrying out DESeq2 analysis now."
    Rscript --vanilla scripts/DESeq2_tiRNA-pipeline.R "$expFile" "$myPath/$outDir/Results/Data/Intermediate-files/"
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
Rscript scripts/Distribution-score.R $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.stddev $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2
### Sort the output but not the header
cat $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.different-distributions.txt | awk 'NR<2{print $0;next}{print $0| "sort -k10,10nr"}' > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.different-distributions.sorted.txt
cat $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.all-features.txt | awk 'NR<2{print $0;next}{print $0| "sort -k10,10nr"}' > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/tiRNA.cond1-vs-cond2.all-features.sorted.txt

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
Rscript scripts/Distribution-score.R $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2.stddev $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2
### Sort the output but not the header
cat $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2.different-distributions.txt | awk 'NR<2{print $0;next}{print $0| "sort -k10,10nr"}' > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2.different-distributions.sorted.txt
cat $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2.all-features.txt | awk 'NR<2{print $0;next}{print $0| "sort -k10,10nr"}' > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/snomiRNA.cond1-vs-cond2.all-features.sorted.txt


### Test whether certain features are cleaved in one condition vs the other



### Get mean and standard deviation (redundant step)
Rscript scripts/Mean_Stdev.R $myPath/$outDir/Results/Data/Intermediate-files/tiRNA.condition1_concatenated.depth $myPath/$outDir/Results/Data/Intermediate-files/condition1_concatenated_mean_stdev.tiRNA.depth &
Rscript scripts/Mean_Stdev.R $myPath/$outDir/Results/Data/Intermediate-files/tiRNA.condition2_concatenated.depth $myPath/$outDir/Results/Data/Intermediate-files/condition2_concatenated_mean_stdev.tiRNA.depth &
Rscript scripts/Mean_Stdev.R $myPath/$outDir/Results/Data/Intermediate-files/snomiRNA.condition1_concatenated.depth $myPath/$outDir/Results/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth &
Rscript scripts/Mean_Stdev.R $myPath/$outDir/Results/Data/Intermediate-files/snomiRNA.condition2_concatenated.depth $myPath/$outDir/Results/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth &
wait

### Get names of differentially expressed features
cat $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/*regulated.csv | grep -v ^,base | awk -F ',' '{print $1}' > $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt
### Get names of features with highest distribution scores 
cat $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/*different-distributions.sorted.txt | grep -v ^feat | awk '{print $1}' > $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/High-distribution-scores_feature-names.txt 
### Get unique set of differentially expressed features and features with high coefficient of variation
cat $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/DESeq2/DEGs_names-only.txt $myPath/$outDir/Results/Data/Intermediate-files/Distribution-score/High-distribution-scores_feature-names.txt | sort | uniq > $myPath/$outDir/Results/Data/DEGs-and-High-distribution-scores_names-only.txt

### If there are differentially expressed/high distribution features, plot these:
if [[ $(wc -l < $myPath/$outDir/Results/Data/DEGs-and-High-distribution-scores_names-only.txt) -ge 1 ]]; then
	### Plot DEGs arg1 and 2 are inputs, arg 3 is list of differentially expressed genes, arg 4 is output pdf, 
	### arg 5 is mean coverage cutoff (plot features with coverage above this), arg 5 is GTF file for snomiRNAs (arg 5 is not given to tiRNA data)
	Rscript scripts/Bedgraph_plotter_DEGs-v3.R $myPath/$outDir/Results/Data/Intermediate-files/condition1_concatenated_mean_stdev.tiRNA.depth $myPath/$outDir/Results/Data/Intermediate-files/condition2_concatenated_mean_stdev.tiRNA.depth $myPath/$outDir/Results/Data/DEGs-and-High-distribution-scores_names-only.txt $myPath/$outDir/Results/Plots/${condition1}_vs_${condition2}_tiRNAs_DEGs.pdf 0 &
	Rscript scripts/Bedgraph_plotter_DEGs-v3.R $myPath/$outDir/Results/Data/Intermediate-files/condition1_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Results/Data/Intermediate-files/condition2_concatenated_mean_stdev.snomiRNA.depth $myPath/$outDir/Results/Data/DEGs-and-High-distribution-scores_names-only.txt $myPath/$outDir/Results/Plots/${condition1}_vs_${condition2}_snomiRNAs_DEGs.pdf 0 DBs/hg19-snomiRNA_cdhit.gtf &
else
	echo "There were no differentially expressed features. No plots were generated." >> $myPath/$outDir/Results/Plots/${condition1}_vs_${condition2}_no-features-to-plot.txt
fi
wait

### Move DESeq results to Data directory
mv $myPath/$outDir/Results/Data/Intermediate-files/DE_Results/ $myPath/$outDir/Results/Data/

#rm -rf Data

echo "Finished at $(date)"

