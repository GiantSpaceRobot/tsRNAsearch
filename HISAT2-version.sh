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
Options:

	-h	Print the usage and options information
	-s	Single-end file for analysis
	-1	First paired-end file
	-2	Second paired-end file
	-o	Output directory for the results and log files
	-T	Use Tophat2 to carry out alignment steps (default: HISAT2)
	-p	Number of CPUs to use {default is to calculate the number of processors and use 75%}

	Input file format should be FASTQ (.fq/.fastq) or gzipped FASTQ (.gz)
	" 1>&2; }

if [ $# -eq 0 ]; then
    echo "No arguments provided"
    asciiArt
	usage
	info
	exit 1
fi

while getopts ":hTp:s:1:2:o:" o; do
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
		p)
			CPUs="$OPTARG"
			;;
		T)
			oldAligner="yes"
			;;
		*)
            echo "Error in input parameters!"
			usage
			info
			exit 1
            ;;
    esac
done

######### Make sure that the absolute paths of files and directories have been specified
### If the pathname specified by $outDir does not begin with a slash, quit (we need full path name)
if [[ ! $outDir = /* ]]; then
    echo "Error: File paths must absolute. Please specify the full path for the output directory."
    exit 1
fi
### If the pathname specified by $expFile does not begin with a slash, quit (we need full path name)
if [ "$singleFile" ]; then
    if [[ ! $singleFile = /* ]]; then
        echo "Error: File paths must absolute. Please specify the full path for the FASTQ input file."
        exit 1
    fi
fi
### If the pathname specified by $expFile does not begin with a slash, quit (we need full path name)
if [ "$file1" ]; then
    if [[ ! $file1 = /* ]]; then
        echo "Error: File paths must absolute. Please specify the full path for the FASTQ input files."
        exit 1
    fi
fi
### If the pathname specified by $expFile does not begin with a slash, quit (we need full path name)
if [ "$file2" ]; then
    if [[ ! $file2 = /* ]]; then
        echo "Error: File paths must absolute. Please specify the full path for the FASTQ input files."
        exit 1
    fi
fi
#########

echo "Started at $(date)" # Print pipeline start-time

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

function bam_to_plots () {  ### Steps for plotting regions with high variation in coverage

	### Output coverage of all features we are interested in (e.g. tRNAs)
	bedtools genomecov -d -split -ibam $1/accepted_hits.bam > $1/accepted_hits.genomecov
	### If we are working with tRNAs, collapse all tRNAs based on same isoacceptor
	if [ $3 = "tiRNA" ]; then
		python scripts/Bedgraph_collapse-tRNAs.py $1/accepted_hits.genomecov $1/accepted_hits_collapsed.genomecov
		mv $1/accepted_hits.genomecov $1/accepted_hits_original.genomecov 
		cp $1/accepted_hits_collapsed.genomecov $1/accepted_hits.genomecov	
	fi
	### Sort by tRNA isoacceptor and nucleotide position
	sort -k1,1 -k2,2n $1/accepted_hits.genomecov > $1/accepted_hits_sorted.genomecov
	### Plot the coverage of all tRNA isoacceptors
	Rscript scripts/Bedgraph_plotter-v3.R $1/accepted_hits_sorted.genomecov $1/$2_$3_Coverage-plots.pdf
	### Output the mean, standard deviation and coefficient of variance of each isoacceptor
	python scripts/Bedgraph-analyser.py $1/accepted_hits_sorted.genomecov $1/accepted_hits_sorted.tsv
	### Gather all tRNA isoacceptors with at least a mean coverage of 10
	awk '$2>10' $1/accepted_hits_sorted.tsv > $1/accepted_hits_sorted_mean-std.tsv
	### Sort the remaining tRNA isoacceptors by coef. of variance
	sort -k4,4nr $1/accepted_hits_sorted_mean-std.tsv > $1/$2_$3_accepted_hits_sorted_mean-std_sorted.tsv
	### Move finalised data for further analysis
	mkdir -p $outDir/Data_and_Plots
	cp $1/$2_$3_Coverage-plots.pdf $outDir/Data_and_Plots/$2_$3_Coverage-plots.pdf
	cp $1/accepted_hits_sorted.genomecov $outDir/Data_and_Plots/$2_$3.genomecov
	cp $1/$2_$3_accepted_hits_sorted_mean-std_sorted.tsv $outDir/Data_and_Plots/$2_$3_genomecov.stats
}


# Estimate CPUs to use if a number has no been provided
if [ -z "$CPUs" ]; then 
	# Determine max number of CPUs available
	maxCPUs=$(grep -c ^processor /proc/cpuinfo)
	# If the variable "CPUs" is not an integer, default to use a value of 1
	re='^[0-9]+$'
	if ! [[ $maxCPUs =~ $re ]] ; then
		echo "Error: Variable 'maxCPUs' is not an integer" >&2
		maxCPUs=1
	fi
	CPUs=$(echo "$maxCPUs * 0.75" | bc | awk '{print int($1+0.5)}') # Use 75% of max CPU number
fi

# Determine if paired-end or single-end files were provided as input
if [ -z "$singleFile" ]; then # If the singleEnd variable is empty
	pairedEnd="True"
else
	pairedEnd="False"
fi

# Check if the output directory exists. If not, create it
string_padder "Creating directory structure"
mkdir -p $outDir
#mkdir $outDir/checkpoints
mkdir -p $outDir/trim_galore_output
mkdir -p $outDir/FastQC
mkdir -p $outDir/tRNA-alignment
mkdir -p $outDir/snomiRNA-alignment
mkdir -p $outDir/mRNA-ncRNA-alignment
mkdir -p $outDir/HTSeq-count-output

#touch $outDir/checkpoints/checkpoint-1.flag

# Run Trim_Galore (paired-end or single-end)
string_padder "Trimming reads using Trim Galore"
if [[ $pairedEnd = "True" ]]; then
	
	echo "
	Input: paired-end read files
	"
	
	string_padder "This version of the pipeline is under construction. Please treat your paired-end files as single-end for now."
	exit 1

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
	string_padder "Trimming complete. Moving on to FastQC analysis"

	# Run FastQC on newly trimmed files
	fastqc -o $outDir/FastQC/ -f fastq $outDir/trim_galore_output/$trimmedFile1 $outDir/trim_galore_output/$trimmedFile2
	string_padder "Finished running FastQC. Moving onto alignment to tRNA database"

	# Align trimmed reads to tRNA database using HISAT2/Tophat2
	if [[ $oldAligner = "yes" ]]; then
		tophat2 -p $CPUs -x 1 -o $outDir/tRNA-alignment/ DBs/bowtie2_index/hg19-wholetRNA-CCA $outDir/trim_galore_output/$trimmedFile1 $outDir/trim_galore_output/$trimmedFile2
		# Tophat2 using HG38 genome below:
		#tophat2 -p $CPUs -x 1 -o $outDir/tRNA-alignment/ DBs/bowtie2_index/hg38-tRNAs_CCA $outDir/trim_galore_output/$trimmedFile1 $outDir/trim_galore_output/$trimmedFile2
		bedtools bamtofastq -i $outDir/tRNA-alignment/unmapped.bam -fq $outDir/tRNA-alignment/unmapped_1.fastq -fq2 $outDir/tRNA-alignment/unmapped_2.fastq
	else
		echo $trimmedFile1
		hisat2 -p $CPUs -x DBs/hisat2_index/hg38-tRNAs_CCA -1 $outDir/trim_galore_output/$trimmedFile1 -2 $outDir/trim_galore_output/$trimmedFile2 -S $outDir/tRNA-alignment/aligned_tRNAdb.sam
        string_padder "This version of the pipeline is not finished yet. Please use the '-T' parameter to use the Tophat2 version"
	fi

elif [[ $pairedEnd = "False" ]]; then

	echo "
	Input: single-end read file
	"
	
	singleFile_base=${singleFile##*/}    # Remove pathname
	singleFile_basename="$( cut -d '.' -f 1 <<< "$singleFile_base" )" # Get filename before first occurence of .
	#suffix="$( cut -d '.' -f 2- <<< "$singleFile_base" )" # Get full file suffix/entension
	if [[ $singleFile == *".gz"* ]]; then
  		suffix="fq.gz"
	else
		suffix="fq"
	fi
	#suffix="fq"
	printf -v trimmedFile "%s_trimmed.%s" "$singleFile_basename" "$suffix"
	printf -v fastqcFile "%s_trimmed_fastqc.html" "$singleFile_basename"
	
	# Run Trim_Galore on single read file
	if [ ! -f $outDir/trim_galore_output/$trimmedFile ]; then
		trim_galore --stringency 10 --length 15 -o $outDir/trim_galore_output/ $singleFile
		string_padder "Trimming complete. Moving on to FastQC analysis"
	else
		string_padder "Found trimmed file. Skipping this step"
	fi
	
	# Run FastQC on newly trimmed file
	if [ ! -f $outDir/FastQC/$fastqcFile ]; then
		# Put fastqcFile definition here?
		fastqc -o $outDir/FastQC/ -f fastq $outDir/trim_galore_output/$trimmedFile
		string_padder "Finished running FastQC. Moving onto alignment to tRNA database"
	else
		string_padder "Found FastQC file. Skipping this step."
	fi

	#touch $outDir/checkpoints/checkpoint-2.flag

	# Align trimmed reads to tRNA database using HISAT2/Tophat2
	
	if [[ $oldAligner = "yes" ]]; then

		# Tophat2 using HG38 genome below
		#tophat2 -p $CPUs -x 1 -o $outDir/tRNA-alignment/ DBs/bowtie2_index/hg38-tRNAs_CCA $outDir/trim_galore_output/$trimmedFile

		if [ ! -f $outDir/tRNA-alignment/align_summary.txt ]; then
			if [ ! -d DBs/tRNA_DB ]; then
				string_padder "Creating one-time tRNA transcriptome"
				tophat2 -p $CPUs -G DBs/hg19-wholetRNA-CCA.gtf --transcriptome-index=DBs/tRNA_DB/known DBs/bowtie2_index/hg19-wholetRNA-CCA
			fi
			string_padder "Running tRNA alignment step..."
			tophat2 -p $CPUs -x 1 -T -N 1 -o $outDir/tRNA-alignment --transcriptome-index=DBs/tRNA_DB/known DBs/bowtie2_index/hg19-wholetRNA-CCA $outDir/trim_galore_output/$trimmedFile
			if [ -f $outDir/tRNA-alignment/unmapped.bam ]; then  #If tophat2 successfully mapped reads, convert the unmapped to FASTQ
				bedtools bamtofastq -i $outDir/tRNA-alignment/unmapped.bam -fq $outDir/tRNA-alignment/$trimmedFile	
				samtools index $outDir/tRNA-alignment/accepted_hits.bam
			else
				echo "
				tRNA alignment output not found. Reads likely did not map to tRNA reference. 
				Using trimmed reads from Trim_Galore output.
				"
				mv $outDir/trim_galore_output/$trimmedFile $outDir/tRNA-alignment/$trimmedFile
			fi
		else
	        string_padder "Found tRNA alignment file. Skipping this step."
		fi

		if [ ! -f $outDir/snomiRNA-alignment/$trimmedFile ]; then # If this file was not generated, try and align the unmapped reads from the tRNA alignment
			if [ ! -d DBs/snomiRNA_DB ]; then
				string_padder "Creating one-time sno/miRNA transcriptome"
				tophat2 -p $CPUs -G DBs/hg19-snomiRNA.gtf --transcriptome-index=DBs/snomiRNA_DB/known DBs/bowtie2_index/Homo_sapiens.GRCh37.dna.primary_assembly
			fi
			string_padder "Runnng sno/miRNA alignment step..."
			tophat2 -p $CPUs -x 1 -T -N 1 -o $outDir/snomiRNA-alignment/ --transcriptome-index=DBs/snomiRNA_DB/known DBs/bowtie2_index/Homo_sapiens.GRCh37.dna.primary_assembly $outDir/tRNA-alignment/$trimmedFile		
		fi
		if [ ! -f $outDir/snomiRNA-alignment/align_summary.txt ]; then # If this file was still not generated, use the unmapped FASTQ from the first alignment output
			echo "
				snomiRNA alignment output not found. Reads likely did not map to sno/miRNA reference. 
				Using unmapped.fastq file from tRNA alignment output.
				"
			mv $outDir/tRNA-alignment/$trimmedFile $outDir/snomiRNA-alignment/$trimmedFile
		else
			samtools index $outDir/snomiRNA-alignment/accepted_hits.bam
			bedtools bamtofastq -i $outDir/snomiRNA-alignment/unmapped.bam -fq $outDir/snomiRNA-alignment/$trimmedFile
		fi
		
		#touch $outDir/checkpoints/checkpoint-3.flag
		if [ ! -f $outDir/mRNA-ncRNA-alignment/align_summary.txt ]; then # If this file was not generated, try and align the unmapped reads from the tRNA alignment	
			if [ ! -d DBs/mRNA-ncRNA_DB ]; then
				string_padder "Creating one-time sno/miRNA transcriptome"
				tophat2 -p $CPUs -G DBs/hg19-mRNA-ncRNA.gtf --transcriptome-index=DBs/mRNA-ncRNA_DB/known DBs/bowtie2_index/Homo_sapiens.GRCh37.dna.primary_assembly
			fi
			string_padder "Running mRNA/ncRNA alignment step..."
			tophat2 -p $CPUs -T -o $outDir/mRNA-ncRNA-alignment/ --transcriptome-index=DBs/mRNA-ncRNA_DB/known DBs/bowtie2_index/Homo_sapiens.GRCh37.dna.primary_assembly $outDir/snomiRNA-alignment/$trimmedFile		
			if [ -f $outDir/mRNA-ncRNA-alignment/unmapped.bam ]; then  #If tophat2 successfully mapped reads, convert the unmapped to FASTQ
				samtools index $outDir/mRNA-ncRNA-alignment/accepted_hits.bam
				bedtools bamtofastq -i $outDir/mRNA-ncRNA-alignment/unmapped.bam -fq $outDir/mRNA-ncRNA-alignment/$trimmedFile
			else
				echo "
				mRNA/ncRNA alignment output not found. Reads likely did not map to mRNA/ncRNA reference. 
				"
				mv $outDir/snomiRNA-alignment/$trimmedFile $outDir/mRNA-ncRNA-alignment/$trimmedFile
				
				## What works here:
				#string_padder "Look at this:"
				#echo $(zcat $outDir/mRNA-ncRNA-alignment/$trimmedFile | wc -l)/4 | bc
				#unmappedReadCount="$(zcat $outDir/mRNA-ncRNA-alignment/$trimmedFile | wc -l)/4|bc"
				##
				
				string_padder "$unmappedReadCount"
        	fi
		else
			string_padder "Found mRNA/ncRNA alignment file. Skipping this step"
		fi		



	else
		#hisat2 -p $CPUs -x DBs/hisat2_index/hg38-tRNAs_CCA -U $outDir/trim_galore_output/$trimmedFile -S $outDir/tRNA-alignment/aligned_tRNAdb.sam
		#string_padder "This version of the pipeline is not finished yet. Please use the '-T' parameter to use the Tophat2 version"
		if [ ! -f $outDir/tRNA-alignment/align_summary.txt ]; then
			#if [ ! -d DBs/tRNA_DB ]; then
			#		string_padder "Creating one-time tRNA transcriptome"
			#	tophat2 -p $CPUs -G DBs/hg19-wholetRNA-CCA.gtf --transcriptome-index=DBs/tRNA_DB/known DBs/bowtie2_index/hg19-wholetRNA-CCA
			#fi
			string_padder "Running tRNA alignment step..."
			#tophat2 -p $CPUs -x 1 -T -N 1 -o $outDir/tRNA-alignment --transcriptome-index=DBs/tRNA_DB/known DBs/bowtie2_index/hg19-wholetRNA-CCA $outDir/trim_galore_output/$trimmedFile
			hisat2 -p $CPUs -x DBs/hisat2_index/hg38-tRNAs_CCA -U $outDir/trim_galore_output/$trimmedFile -S $outDir/tRNA-alignment/aligned_tRNAdb.sam --summary-file $outDir/tRNA-alignment/align_summary.txt --un $outDir/tRNA-alignment/unmapped.fastq
			samtools view -bS $outDir/tRNA-alignment/aligned_tRNAdb.sam > $outDir/tRNA-alignment/accepted_hits.bam
		   	echo "Quitting now"	
			exit 1
			if [ -f $outDir/tRNA-alignment/unmapped.bam ]; then  #If tophat2 successfully mapped reads, convert the unmapped to FASTQ
				bedtools bamtofastq -i $outDir/tRNA-alignment/unmapped.bam -fq $outDir/tRNA-alignment/$trimmedFile	
				samtools index $outDir/tRNA-alignment/accepted_hits.bam
			else
				echo "
				tRNA alignment output not found. Reads likely did not map to tRNA reference. 
				Using trimmed reads from Trim_Galore output.
				"
				mv $outDir/trim_galore_output/$trimmedFile $outDir/tRNA-alignment/$trimmedFile
			fi
		else
	        string_padder "Found tRNA alignment file. Skipping this step."
		fi

		if [ ! -f $outDir/snomiRNA-alignment/$trimmedFile ]; then # If this file was not generated, try and align the unmapped reads from the tRNA alignment
			if [ ! -d DBs/snomiRNA_DB ]; then
				string_padder "Creating one-time sno/miRNA transcriptome"
				tophat2 -p $CPUs -G DBs/hg19-snomiRNA.gtf --transcriptome-index=DBs/snomiRNA_DB/known DBs/bowtie2_index/Homo_sapiens.GRCh37.dna.primary_assembly
			fi
			string_padder "Runnng sno/miRNA alignment step..."
			tophat2 -p $CPUs -x 1 -T -N 1 -o $outDir/snomiRNA-alignment/ --transcriptome-index=DBs/snomiRNA_DB/known DBs/bowtie2_index/Homo_sapiens.GRCh37.dna.primary_assembly $outDir/tRNA-alignment/$trimmedFile		
		fi
		if [ ! -f $outDir/snomiRNA-alignment/align_summary.txt ]; then # If this file was still not generated, use the unmapped FASTQ from the first alignment output
			echo "
				snomiRNA alignment output not found. Reads likely did not map to sno/miRNA reference. 
				Using unmapped.fastq file from tRNA alignment output.
				"
			mv $outDir/tRNA-alignment/$trimmedFile $outDir/snomiRNA-alignment/$trimmedFile
		else
			samtools index $outDir/snomiRNA-alignment/accepted_hits.bam
			bedtools bamtofastq -i $outDir/snomiRNA-alignment/unmapped.bam -fq $outDir/snomiRNA-alignment/$trimmedFile
		fi
		
		#touch $outDir/checkpoints/checkpoint-3.flag
		if [ ! -f $outDir/mRNA-ncRNA-alignment/align_summary.txt ]; then # If this file was not generated, try and align the unmapped reads from the tRNA alignment	
			if [ ! -d DBs/mRNA-ncRNA_DB ]; then
				string_padder "Creating one-time sno/miRNA transcriptome"
				tophat2 -p $CPUs -G DBs/hg19-mRNA-ncRNA.gtf --transcriptome-index=DBs/mRNA-ncRNA_DB/known DBs/bowtie2_index/Homo_sapiens.GRCh37.dna.primary_assembly
			fi
			string_padder "Running mRNA/ncRNA alignment step..."
			tophat2 -p $CPUs -T -o $outDir/mRNA-ncRNA-alignment/ --transcriptome-index=DBs/mRNA-ncRNA_DB/known DBs/bowtie2_index/Homo_sapiens.GRCh37.dna.primary_assembly $outDir/snomiRNA-alignment/$trimmedFile		
			if [ -f $outDir/mRNA-ncRNA-alignment/unmapped.bam ]; then  #If tophat2 successfully mapped reads, convert the unmapped to FASTQ
				samtools index $outDir/mRNA-ncRNA-alignment/accepted_hits.bam
				bedtools bamtofastq -i $outDir/mRNA-ncRNA-alignment/unmapped.bam -fq $outDir/mRNA-ncRNA-alignment/$trimmedFile
			else
				echo "
				mRNA/ncRNA alignment output not found. Reads likely did not map to mRNA/ncRNA reference. 
				"
				mv $outDir/snomiRNA-alignment/$trimmedFile $outDir/mRNA-ncRNA-alignment/$trimmedFile
				
				## What works here:
				#string_padder "Look at this:"
				#echo $(zcat $outDir/mRNA-ncRNA-alignment/$trimmedFile | wc -l)/4 | bc
				#unmappedReadCount="$(zcat $outDir/mRNA-ncRNA-alignment/$trimmedFile | wc -l)/4|bc"
				##
				
				string_padder "$unmappedReadCount"
        	fi
		else
			string_padder "Found mRNA/ncRNA alignment file. Skipping this step"
		fi
	
	
	
	
	
	
	
	
	
	fi
fi

#touch $outDir/checkpoints/checkpoint-4.flag

# Produce read counts for the three alignment steps. If one of the alignment steps failed, use an empty htseq-count output file.
string_padder "Alignment steps complete. Moving on to read-counting using HTSeq-count"

# Count for alignment step 1
if [ ! -f $outDir/tRNA-alignment/accepted_hits.bam ]; then
	echo "
No alignment file found for tRNA alignment. Using blank count file instead
"
	cp additional-files/empty_tRNA.count $outDir/HTSeq-count-output/tRNA-alignment.count
else
	echo "
Counting tRNA alignment reads
"
	htseq-count -f bam $outDir/tRNA-alignment/accepted_hits.bam DBs/hg19-wholetRNA-CCA.gtf > $outDir/HTSeq-count-output/tRNA-alignment.count
	bam_to_plots $outDir/tRNA-alignment $singleFile_basename tiRNA
fi

# Count for alignment step 2
if [ ! -f $outDir/snomiRNA-alignment/accepted_hits.bam ]; then
	echo "
No alignment file found for sno/miRNA alignment. Using blank count file instead
"
	cp additional-files/empty_snomiRNA.count $outDir/HTSeq-count-output/snomiRNA-alignment.count
else
	echo "
Counting sno/miRNA alignment reads
"
	htseq-count -f bam $outDir/snomiRNA-alignment/accepted_hits.bam DBs/hg19-snomiRNA.gtf > $outDir/HTSeq-count-output/snomiRNA-alignment.count
    bam_to_plots $outDir/snomiRNA-alignment $singleFile_basename snomiRNA
fi

#touch $outDir/checkpoints/checkpoint-5.flag
# Count for alignment step 3
if [ ! -f $outDir/mRNA-ncRNA-alignment/accepted_hits.bam ]; then
	echo "
No alignment file found for mRNA/ncRNA alignment. Using blank count file instead
"
	cp additional-files/empty_mRNA-ncRNA.count $outDir/HTSeq-count-output/mRNA-ncRNA-alignment.count
else
	echo "
Counting mRNA/ncRNA alignment reads
"
	htseq-count -f bam $outDir/mRNA-ncRNA-alignment/accepted_hits.bam DBs/Homo_sapiens.GRCh37.87.gtf > $outDir/HTSeq-count-output/mRNA-ncRNA-alignment.count
    #bam_to_plots $outDir/mRNA-ncRNA-alignment # The resulting file would be too big. It would be interesting to see the top few genes/ncRNAs and their coverage.
fi

#touch $outDir/checkpoints/checkpoint-6.flag


# Create a log file with the date, time and name of the input file in it's name


