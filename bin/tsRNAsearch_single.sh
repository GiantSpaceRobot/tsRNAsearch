#!/bin/bash

# Usage: tsRNAsearch_single.sh -h
# Author: Paul Donovan
# Email: pauldonovan@rcsi.com
# 19-Oct-2018

asciiArt() { echo '

 _        _____ _   _   ___                          _     
| |      | ___ \ \ | | / _ \                        | |    
| |_ ___ | |_/ /  \| |/ /_\ \___  ___  __ _ _ __ ___| |__  
| __/ __||    /| . ` ||  _  / __|/ _ \/ _` | `__/ __| `_ \ 
| |_\__ \| |\ \| |\  || | | \__ \  __/ (_| | | | (__| | | |
 \__|___/\_| \_\_| \_/\_| |_/___/\___|\__,_|_|  \___|_| |_|
                                                          

	' 1>&1; }
usage() { echo "Usage: $0 -o OutputDirectory/ -f /path/to/SeqFile.fastq.gz
" 1>&2; }
info() { echo "
Options:

	-h	Print the usage and options information
	-s	Analyse data against 'human', 'mouse', or 'rat'? {default: human}
	-f	Single-end file for analysis
	-o	Output directory for the results and log files
	-A	Plot all features? yes/no {default: yes}
	-t	Number of threads to use {default is to calculate the number of processors and use 75%}
	-S	Skip pre-processing of data (i.e. skip FastQC and trim_galore) {default: no}

	Input file format should be FASTQ (.fq/.fastq) or gzipped FASTQ (.gz)
	" 1>&2; }

if [ $# -eq 0 ]; then
    echo "No arguments provided"
    asciiArt
	usage
	info
	exit 1
fi

### Define defaults
species="human"
skip="no"
Plots="yes"
while getopts ":hs:t:f:o:A:S:R:" o; do
    case "${o}" in
		h)
			asciiArt
			usage
			info
			exit
			;;
		s)
			species="$OPTARG"
			;;
		f)
			singleFile="$OPTARG"
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
		R)
			remove="$OPTARG"
			;;
		S)
			skip="$OPTARG"
			;;
		*)
            echo "Error in input parameters!"
			usage
			info
			exit 1
            ;;
    esac
done

echo "Started analysing "$singleFile" on $(date)" # Print pipeline start-time

ncRNADB="DBs/species_index/${species}-ncRNAs/"
tRNAGTF="DBs/${species}_tRNAs_relative_cdhit.gtf"
ncRNA_GTF="DBs/${species}_ncRNAs_relative_cdhit.gtf"
tRNA_introns="additional-files/${species}_tRNA-introns-for-removal.tsv"
empty_tRNAs="additional-files/${species}_empty_tRNA.count"
empty_ncRNAs="additional-files/${species}_empty_ncRNA.count"

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

function SAMcollapse () {
	string_padder "Collapsing SAM file..."
	### How many chunks to create:
	readNumber=$(grep -v ^@ $1 | awk '{print $1}' | uniq | wc -l)
	if (( $readNumber > 100000 )); then
		chunksRaw=$(( readNumber / 10000 )) # 10,000 lines per SAM chunk
		chunks=$( echo $chunksRaw | awk '{print int($1+0.5)}' )
		echo "Over 100,000 unique reads in SAM file. Splitting SAM into $chunks files..."
	elif (( readNumber>=10000 && readNumber<=100000 )); then
		chunksRaw=$(( readNumber / 5000 )) # 5000 lines per SAM chunk
		chunks=$( echo $chunksRaw | awk '{print int($1+0.5)}' )
		echo "Between 10,000 and 100,000 unique reads in SAM file. Splitting SAM into $chunks files..."
	elif (( readNumber>=1000 && readNumber<=10000 )); then
		chunks=10 #
		echo "Between 1,000 and 10,000 unique reads in SAM file. Splitting SAM into $chunks files..."
	else
		chunks=1
		echo "Less than 1,000 unique reads in SAM file. Splitting SAM into $chunks files..."
	fi
	### Define variables
	threads_available_for_chunks=$threads
	if (( $threads_available_for_chunks > $chunks )); then
		# make sure the number of threads is not set higher than no. of files (chunks) after split
		threads_available_for_chunks=$chunks
	fi
	echo "
		Threads: $threads_available_for_chunks
		# Split SAM files: $chunks"
	fileLen=$(< "$1" wc -l)
	division1=$((fileLen/chunks))
	division=$((division1 + 1))
	myFile="tempFile"
	mkdir -p $outDir/tempDir

	### Remove header
	grep ^@ $1 > $outDir/tempDir/myHeader.txt &
	grep -v ^@ $1 > $outDir/tempDir/mySAM.sam &
	wait
	### 
	fileToCollapse=$outDir/tempDir/mySAM.sam
	### Split file
	echo "Splitting SAM..."
	split -l $division $fileToCollapse $outDir/tempDir/splitFile_
	### Gather first and last read from every split file and add to separate file. Remove these reads from the split files.
	echo "Gathering the names of the first and last reads from every SAM chunk..."
	for i in $outDir/tempDir/splitFile_*; do
		base=$(basename $i)
		first=$(awk 'NR==1' $i | awk '{print $1}') 
		echo $first >> $outDir/tempDir/${myFile}_HeadsAndTails.txt
		last=$(awk 'END{print}' $i | awk '{print $1}') 
		echo $last >> $outDir/tempDir/${myFile}_HeadsAndTails.txt
	done
	echo "Gathering unique set of read names from first/last read names..."
	sort $outDir/tempDir/${myFile}_HeadsAndTails.txt | uniq \
		> $outDir/tempDir/${myFile}_HeadsAndTails_uniq.txt #remove duplicates
	sed -i 's/$/\t/' $outDir/tempDir/${myFile}_HeadsAndTails_uniq.txt # Add tab to end of every line to match pattern exactly
	grep -f $outDir/tempDir/${myFile}_HeadsAndTails_uniq.txt $fileToCollapse \
		> $outDir/tempDir/edit_heads-and-tails #grep all patterns from the heads/tails file
	echo "Extracting alignments for first/last reads from all files..."
	for i in $outDir/tempDir/splitFile_*; do
		base=$(basename $i)
		grep -v -f  $outDir/tempDir/${myFile}_HeadsAndTails_uniq.txt $i \
			> $outDir/tempDir/edit_${base}
	done
	### Run SAMcollapse.py. This loop will only run $threads_available_for_chunks processes at once
	COUNTER=1
	if (( $chunks == 1)); then # If there are so few reads in input that only one SAM chunk is present
		chunksDiv=1
	else # This is the most common situation: Very large SAM has been split into many chunks, divide the chunks by 10 to allow printing of progress
		chunksDiv=$((chunks/10))
	fi
	echo "Collapsing every chunk of SAM..."
	for i in $outDir/tempDir/edit_*; 
	do
		base=$(basename $i)
		if (( $chunks > 1)); then
			python2 bin/SAMcollapse.py \
				$i \
				${fileToCollapse}_${base} \
				>> $outDir/tRNA-alignment/collapsed-reads.txt & 
			if (( $COUNTER % $chunksDiv == 0 )); then # If SAM is split into 100 chunks, print progress when the counter reaches 10, 20, 30... 
				echo "Started job $COUNTER of $chunks" # i.e. print progress every 10%
			fi
			numjobs=($(jobs | wc -l))
			COUNTER=$[$COUNTER + 1]
			while (( $numjobs == $threads_available_for_chunks )); do
				numjobs=($(jobs | wc -l))
				sleep 2 #Enter next loop iteration
			done
		else
			python2 bin/SAMcollapse.py \
				$i \
				${fileToCollapse}_${base} \
				>> $outDir/tRNA-alignment/collapsed-reads.txt & 
			echo "Started job $COUNTER of $chunks" 
		fi
	done
	wait
	readsCollapsedSpecies=$(awk '{split($0,a," "); sum += a[1]} END {print sum}' $outDir/tRNA-alignment/collapsed-reads.txt)
	readsCollapsedGroup=$(awk '{split($0,a," "); sum += a[2]} END {print sum}' $outDir/tRNA-alignment/collapsed-reads.txt)
	echo -e "SAM collapse results:\n\t$readsCollapsedSpecies reads collapsed at the tRNA species level (e.g. 2 gene copies of ProCCG)\n\t$readsCollapsedGroup reads collapsed at the tRNA group level (e.g. ProCCG and ProAAG)"

	### Concatenate results
	echo "Gathering reads that were mapped to similar tRNAs..."
	echo -e "tRNA.group\tread.start\tread.end.approx\tread.name" \
		> $outDir/tRNA-alignment/tRNAs-almost-mapped.txt
	cat $outDir/tempDir/*tRNAs-almost-mapped* | sort \
		>> $outDir/tRNA-alignment/tRNAs-almost-mapped.txt
	mkdir $outDir/tempDir/tRNAsAlmostMapped
	mv $outDir/tempDir/*tRNAs-almost-mapped* $outDir/tempDir/tRNAsAlmostMapped/
	echo "Concatenating SAM header with collapsed files..."
	cat $outDir/tempDir/myHeader.txt ${fileToCollapse}*_edit_* \
		> $outDir/Collapsed.sam
	rm -rf $outDir/tempDir/ # Remove temp directory
	echo "Finished collapsing SAM file"
}

function bam_to_plots () {  ### Steps for plotting regions with high variation in coverage
	### Output coverage of all features we are interested in (e.g. tRNAs)
	samtools depth \
		-d 100000000 \
		-aa $1/accepted_hits.bam \
		> $1/accepted_hits.depth   # A lot faster than bedtools genomecov
	cp $1/accepted_hits.depth $1/accepted_hits_raw.depth
	### Normalise by reads per million (RPM)
	python2 bin/Depth-to-Depth_RPM.py \
		$1/accepted_hits_raw.depth \
		$mapped \
		$1/accepted_hits.depth 
	### If we are working with tRNAs, collapse all tRNAs based on same isoacceptor
	if [[ $3 = "tsRNA" ]]; then
		### Remove introns from tRNA counts (as these will interfere with the read counts of collapsed tRNAs)
		if [ -s $tRNA_introns ]; then # If there are introns in the intron annotation file, proceed with removing them
			Rscript bin/Remove-introns.R \
				$1/accepted_hits.depth \
				$tRNA_introns \
				$1/accepted_hits_intron-removed.depth 
		else # If not, just use the input
			cp $1/accepted_hits.depth $1/accepted_hits_intron-removed.depth
		fi
		### Flip the read coverage of tRNAs that are in the minus orientation
		## This step is not required when using GtRNAdb derived FASTA files as these are all in the plus orientation
		#Rscript bin/Coverage-flipper.R \
		#	$1/accepted_hits_intron-removed.depth \
		#	$tRNAGTF \
		#	$1/accepted_hits_flipped.depth
		cp $1/accepted_hits_intron-removed.depth $1/accepted_hits_flipped.depth
		### Collapse tRNAs from the same tRNA species
		python2 bin/Bedgraph_collapse-tRNAs.py \
			$1/accepted_hits_flipped.depth \
			$1/accepted_hits_collapsed.depth
		### Rename input depth file
		mv $1/accepted_hits.depth $1/accepted_hits_original.depth 
		### Copy the collapsed depth file and name it so that the remaining steps below do not have errors
		cp $1/accepted_hits_collapsed.depth $1/accepted_hits.depth
		### Plot tRNA alignment lengths
		grep -v ^@ $1/tsRNAs_aligned.sam > $1/tsRNAs_aligned.noheader.sam 
		Rscript bin/tRNA_Alignment_Length.R \
				$1/tsRNAs_aligned.noheader.sam \
				$1/$2_$3_tRNA-alignment-length &
	fi
	### Sort by feature name and nucleotide position
	sort -k1,1 -k2,2n $1/accepted_hits.depth \
		> $1/accepted_hits_sorted.depth
	### If -A parameter was provided, plot everything
	if [[ $Plots == "yes" ]]; then # Plot everything
		### Plot the coverage of all features (arg 3 is mean coverage in RPM) and
		### Create plot and txt file describing relationship between 5' and 3' regions of feature
		if [[ $3 = "ncRNA" ]]; then
			Rscript bin/Bedgraph_plotter.R \
				$1/accepted_hits_sorted.depth \
				$1/$2_$3_Coverage-plots.pdf \
				0 \
				$ncRNA_GTF
			Rscript bin/Single-replicate-analysis.R \
				$1/accepted_hits_sorted.depth \
				$1/$2_$3_Results \
				$ncRNA_GTF &
		elif [[ $3 == "tsRNA" ]]; then
			Rscript bin/Bedgraph_plotter.R \
				$1/accepted_hits_sorted.depth \
				$1/$2_$3_Coverage-plots.pdf \
				0
			Rscript bin/Single-replicate-analysis.R \
				$1/accepted_hits_sorted.depth \
				$1/$2_$3_Results &	
		else
			Rscript bin/Bedgraph_plotter.R \
				$1/accepted_hits_sorted.depth \
				$1/$2_$3_Coverage-plots.pdf \
				1000 \
				$speciesGTF
			Rscript bin/Single-replicate-analysis.R \
				$1/accepted_hits_sorted.depth \
				$1/$2_$3_Results &
		fi
		cp $1/$2_$3_Coverage-plots.pdf $outDir/Data_and_Plots/$2_$3_Coverage-plots.pdf
	fi
	### Output the mean, standard deviation and coefficient of variance of each ncRNA/gene
	python2 bin/Bedgraph-analyser.py \
		$1/accepted_hits_sorted.depth \
		$1/accepted_hits_sorted.tsv
	### Gather all ncRNAs/genes with a mean coverage greater than 0 (pointless step but the cutoff used to be higher than 0)
	awk '$2>0' $1/accepted_hits_sorted.tsv \
		> $1/accepted_hits_sorted_mean-std.tsv
	### Sort the remaining ncRNAs/genes by coef. of variance
	sort -k4,4nr $1/accepted_hits_sorted_mean-std.tsv \
		> $1/$2_$3_accepted_hits_sorted_mean-std_sorted.tsv
	### Move finalised data for further analysis
	cp $1/accepted_hits_sorted.depth $outDir/Data_and_Plots/$2_$3.depth &
	cp $1/$2_$3_accepted_hits_sorted_mean-std_sorted.tsv $outDir/Data_and_Plots/$2_$3_depth.inf
	echo -e "Feature\tMean\tStandard Deviation\tCoefficient of Variation\n" > \
		$outDir/Data_and_Plots/Header.txt
	cat $outDir/Data_and_Plots/Header.txt $outDir/Data_and_Plots/$2_$3_depth.inf \
		> $outDir/Data_and_Plots/$2_$3_depth.stats
	rm $outDir/Data_and_Plots/$2_$3_depth.stats $outDir/Data_and_Plots/Header.txt
	wait
}

# If -A parameter was not provided, default is to plot everything
if [ ! "$Plots" ]; then
    Plots="yes"
fi

# Estimate threads to use if a number has no been provided
if [ -z "$threads" ]; then 
	# Determine max number of threads available
	maxthreads=$(grep -c ^processor /proc/cpuinfo)
	# If the variable "threads" is not an integer, default to use a value of 1
	re='^[0-9]+$'
	if ! [[ $maxthreads =~ $re ]] ; then
		echo "Error: Variable 'maxthreads' is not an integer" >&2
		maxthreads=1
	fi
	threads=$(echo "$maxthreads * 0.75" | bc | awk '{print int($1+0.5)}') # Use 75% of max CPU number
fi

### Do not exceed 32 threads for featureCounts as per user manual guidelines
if (( $threads > 32 )); then 
	featureCount_threads=32
else
	featureCount_threads=$threads
fi

# Check if the output directory exists. If not, create it
string_padder "Creating directory structure"
mkdir -p $outDir
mkdir -p $outDir/FastQC
mkdir -p $outDir/tRNA-alignment
mkdir -p $outDir/ncRNA-alignment
mkdir -p $outDir/FCount-count-output
mkdir -p $outDir/FCount-to-RPM
mkdir -p $outDir/Data_and_Plots
mkdir -p $outDir/Data_and_Plots/Encoded

singleFile_base=${singleFile##*/}    # Remove pathname
singleFile_basename="$( cut -d '.' -f 1 <<< "$singleFile_base" )" # Get filename before first occurence of .
if [[ $singleFile == *".gz"* ]]; then
	suffix="fq.gz"
	STARparam="--readFilesCommand zcat"
else
	suffix="fq"
	STARparam=""
fi
printf -v trimmedFile "%s_trimmed.%s" "$singleFile_basename" "$suffix"
printf -v myFile "%s.fq" "$singleFile_basename"
printf -v fastqcFile "%s_trimmed_fastqc.html" "$singleFile_basename"

### Skip pre-processing or not:
if [ $skip = "no" ]; then
	# Make directories for pre-processing
	mkdir -p $outDir/pre-processing
	string_padder "Pre-processing reads using trim_galore..."
	bin/trim_galore \
		--stringency 10 \
		--length 15 \
		-o $outDir/pre-processing/ \
		--fastqc_args "--outdir $outDir/FastQC/" \
		$singleFile
	readsForAlignment=$outDir/pre-processing/$trimmedFile
else
	# Skip pre-processing by using the provided RNA-seq file
	readsForAlignment=$singleFile
fi

### Align reads to ncRNA database using STAR
string_padder "Running tRNA/ncRNA alignment step..."

### STAR ###
bin/STAR \
	--runThreadN $threads \
	--genomeDir $ncRNADB \
	--readFilesIn $readsForAlignment \
	--outFileNamePrefix $outDir/tRNA-alignment/ \
	--outSAMattributes AS nM HI NH \
	--outFilterMultimapScoreRange 0 \
	--outReadsUnmapped Fastx \
	--outFilterMatchNmin 15 \
	$STARparam
grep "Number of input reads" $outDir/tRNA-alignment/Log.final.out \
	| awk -F '\t' '{print $2}' \
	| tr -d '\040\011\012\015' \
	> $outDir/Stats.log # the tr command removes all types of spaces
echo " reads; of these:" >> $outDir/Stats.log
## Collapse SAM files:
SAMcollapse $outDir/tRNA-alignment/Aligned.out.sam #Collapse reads aligned to the same tRNA species 
mv $outDir/Collapsed.sam $outDir/tRNA-alignment/aligned_tRNAdb.sam # match hisat2 naming convention
mv $outDir/tRNA-alignment/Unmapped.out.mate1 $outDir/tRNA-alignment/unmapped.fastq
## Or don't collapse SAM files (comment out three lines above and uncomment line below)
#mv $outDir/tRNA-alignment/Aligned.out.sam $outDir/tRNA-alignment/aligned_tRNAdb.sam
#### STAR ###

if [ -f $outDir/tRNA-alignment/aligned_tRNAdb.sam ]; then  #If STAR successfully mapped reads, convert to bam and index
	### Split the resulting SAM file into reads aligned to tRNAs and ncRNAs
	grep ^@ $outDir/tRNA-alignment/aligned_tRNAdb.sam > $outDir/tRNA-alignment/SamHeader.sam &
	grep ENS $outDir/tRNA-alignment/aligned_tRNAdb.sam > $outDir/tRNA-alignment/ncRNAs.sam &
	grep -v ENS $outDir/tRNA-alignment/aligned_tRNAdb.sam > $outDir/tRNA-alignment/tsRNAs_aligned.sam &
	wait
	cat $outDir/tRNA-alignment/SamHeader.sam $outDir/tRNA-alignment/ncRNAs.sam \
		> $outDir/tRNA-alignment/ncRNAs_aligned.sam
	wait
	samtools view -bS $outDir/tRNA-alignment/tsRNAs_aligned.sam \
		| samtools sort \
		> $outDir/tRNA-alignment/accepted_hits.bam
	samtools index $outDir/tRNA-alignment/accepted_hits.bam &
	mv $outDir/tRNA-alignment/unmapped.fastq $outDir/tRNA-alignment/$myFile
	### Move ncRNA BAM to directory
	samtools view -bS $outDir/tRNA-alignment/ncRNAs.sam \
		| samtools sort \
		> $outDir/ncRNA-alignment/accepted_hits.bam
	samtools index $outDir/ncRNA-alignment/accepted_hits.bam &
	rm \
		$outDir/tRNA-alignment/SamHeader.sam \
		$outDir/tRNA-alignment/ncRNAs.sam \
		$outDir/tRNA-alignment/aligned_tRNAdb.sam &	
else
	echo "
	Alignment output not found. Reads likely did not map to ncRNA database. 
	Using trimmed reads from Trim_Galore output.
	"
	cp $outDir/trim_galore_output/$trimmedFile $outDir/tRNA-alignment/$trimmedFile
fi

# Produce read counts for the three alignment steps. If one of the alignment steps failed, use an empty htseq-count output file.
string_padder "Alignment steps complete. Moving on to read-counting using FCount-count"

# Count for alignment step 2
if [ ! -f $outDir/ncRNA-alignment/accepted_hits.bam ]; then
	echo "
No alignment file found for ncRNA alignment. Using blank count file instead
"
	cp $empty_ncRNAs $outDir/FCount-count-output/ncRNA-alignment.count &
else
	echo "
Counting ncRNA alignment reads
"
	bin/featureCounts \
		-T $featureCount_threads \
		-a $ncRNA_GTF \
		-o $outDir/FCount-count-output/ncRNA-alignment.fcount \
		$outDir/ncRNA-alignment/accepted_hits.bam
	grep -v featureCounts $outDir/FCount-count-output/ncRNA-alignment.fcount \
		| grep -v ^Geneid \
		| awk -v OFS='\t' '{print $1, $7}' \
		> $outDir/FCount-count-output/ncRNA-alignment.count
fi

# Count for alignment step 1
if [ ! -f $outDir/tRNA-alignment/accepted_hits.bam ]; then
	echo "
No alignment file found for tRNA alignment. Using blank count file instead
"
	cp $empty_tRNAs $outDir/FCount-count-output/tRNA-alignment.count &
else
	echo "
Counting tRNA alignment reads
"
	bin/featureCounts \
		-T $featureCount_threads \
		-a $tRNAGTF \
		-o $outDir/FCount-count-output/tRNA-alignment.fcount \
		$outDir/tRNA-alignment/accepted_hits.bam
	grep -v featureCounts $outDir/FCount-count-output/tRNA-alignment.fcount \
		| grep -v ^Geneid \
		| awk -v OFS='\t' '{print $1, $7}' \
		> $outDir/FCount-count-output/tRNA-alignment.count
fi

wait

### Get total reads mapped
string_padder "Get total number of reads mapped"
cat $outDir/FCount-count-output/*.count | grep -v ^__ | sort -k1,1 \
	> $outDir/FCount-count-output/$singleFile_basename.all-features.count 
sed -i '1s/^/Features\t'"$singleFile_basename"'\n/' $outDir/FCount-count-output/$singleFile_basename.all-features.count # Add column headers
mapped=$(awk '{sum+=$2} END{print sum;}' $outDir/FCount-count-output/$singleFile_basename.all-features.count)
echo "Reads mapped: $mapped" >> $outDir/Stats.log
echo "Reads mapped: $mapped"
wait

### Plot everything
string_padder "Generate depth files and plot features"
bam_to_plots $outDir/tRNA-alignment $singleFile_basename tsRNA &
bam_to_plots $outDir/ncRNA-alignment $singleFile_basename ncRNA &

### Process multi-mapping tRNAs
python2 bin/Leftovers-to-Bedgraph.py \
	$outDir/tRNA-alignment/tRNAs-almost-mapped.txt \
	additional-files/${species}_tRNA-lengths.txt \
	$outDir/tRNA-alignment/tRNAs-almost-mapped.depth
python2 bin/Depth-to-Depth_RPM.py \
	$outDir/tRNA-alignment/tRNAs-almost-mapped.depth \
	$mapped \
	$outDir/tRNA-alignment/$singleFile_basename.tRNAs-almost-mapped_RPM.depth
sort -k1,1 -k2,2n $outDir/tRNA-alignment/$singleFile_basename.tRNAs-almost-mapped_RPM.depth \
	> $outDir/tRNA-alignment/$singleFile_basename.tRNAs-almost-mapped_RPM.sorted.depth
mv $outDir/tRNA-alignment/$singleFile_basename.tRNAs-almost-mapped_RPM.sorted.depth \
	$outDir/tRNA-alignment/$singleFile_basename.tRNAs-almost-mapped_RPM.depth

### Get RPM-normalised FCount count data
string_padder "Get RPM-normalised read-counts"
python2 bin/FCount-to-RPM.py \
	$outDir/FCount-count-output/$singleFile_basename.all-features.count \
	$mapped \
	$outDir/FCount-to-RPM/$singleFile_basename.all-features &
python2 bin/FCount-to-RPM.py \
	$outDir/FCount-count-output/tRNA-alignment.count \
	$mapped \
	$outDir/FCount-to-RPM/tRNA-alignment & 
python2 bin/FCount-to-RPM.py \
	$outDir/FCount-count-output/ncRNA-alignment.count \
	$mapped \
	$outDir/FCount-to-RPM/ncRNA-alignment &
wait
sleep 5  # Make sure everything is finished running

### Collapse count file
string_padder "Collapsing count files..."
python2 bin/CollapseCountfile.py \
	$outDir/FCount-count-output/$singleFile_basename.all-features.count \
	$outDir/FCount-count-output/$singleFile_basename.collapsed.all-features.count  # For DESeq2
python2 bin/CollapseCountfile.py \
	$outDir/FCount-to-RPM/$singleFile_basename.all-features.rpm.count \
	$outDir/FCount-to-RPM/$singleFile_basename.collapsed.all-features.rpm.count # For Cleavage + Distribution algorithms
wait

### Move results to Data_and_Plots
cp $outDir/FCount-to-RPM/$singleFile_basename.collapsed.all-features.rpm.count $outDir/Data_and_Plots/
### Move tRNA alignment length plot
cp $outDir/tRNA-alignment/${singleFile_basename}_tsRNA_tRNA-alignment-length.pdf $outDir/Data_and_Plots/
if [[ $Plots == "yes" ]]; then
	### If extra plotting parameter (-A) was selected, copy these files 
	Rscript bin/Bedgraph_plotter.R \
		$outDir/tRNA-alignment/$singleFile_basename.tRNAs-almost-mapped_RPM.depth \
		$outDir/tRNA-alignment/Multi-mappers_tsRNAs_Coverage-plots.pdf \
		0
	cp $outDir/tRNA-alignment/Multi-mappers_tsRNAs_Coverage-plots.pdf $outDir/Data_and_Plots/
	cp $outDir/tRNA-alignment/*Results.* $outDir/Data_and_Plots/
	cp $outDir/ncRNA-alignment/*Results.* $outDir/Data_and_Plots/
fi

### Determine tsRNA types for each tRNA profile
grep -v ^feature $outDir/tRNA-alignment/tRNAs-almost-mapped.depth \
	> $outDir/tRNA-alignment/tRNAs-almost-mapped_no-header.depth
cat $outDir/Data_and_Plots/${singleFile_basename}_tsRNA.depth \
	$outDir/tRNA-alignment/tRNAs-almost-mapped_no-header.depth \
	> $outDir/tRNA-alignment/All-tRNAs.depth
Rscript bin/tsRNA-type-classification.R \
	$outDir/tRNA-alignment/All-tRNAs.depth \
	$outDir/Data_and_Plots/tsRNAs-classified-by-type.txt 
awk '{print $1"\t"$2}' $outDir/Data_and_Plots/tsRNAs-classified-by-type.txt \
	> $outDir/tRNA-alignment/tsRNAs-classified-by-type_clean.txt
sed -e 's/^/<br \/>/' $outDir/tRNA-alignment/tsRNAs-classified-by-type_clean.txt \
	> $outDir/tRNA-alignment/tsRNAs-classified-by-type_clean_HTML.txt
tsRNAtype=$(cat $outDir/tRNA-alignment/tsRNAs-classified-by-type_clean_HTML.txt)
sed -e 's/^/* /' $outDir/tRNA-alignment/tsRNAs-classified-by-type_clean.txt \
	> $outDir/tRNA-alignment/tsRNAs-classified-by-type_clean_RmdHTML.txt

### If -R == yes, remove intermediate files
if [[ $remove == "yes" ]]; then 
	find $outDir/ -name '*bam' -delete   # Remove BAMs
	find $outDir/ -name '*sam' -delete   # SAMs
	find $outDir/ -name '*tempFile' -delete   # tempFiles
	find $outDir/ -name '*.gz' -delete   # .gz files
	find $outDir/ -name '*fastq' -delete # .fastq files
	find $outDir/ -name '*fq' -delete    # .fq files
	find $outDir/ -name '*accepted_hits*depth' -delete    # .depth files
fi

### Encode figures using base64
for i in $outDir/Data_and_Plots/*; do # Encode using base64
	my_file=$(basename $i)
	if [[ $my_file =~ ".pdf" ]]; then # Encode PDFs
		echo -n "<iframe src=\"data:application/pdf;base64," > $outDir/Data_and_Plots/Encoded/$my_file.txt
		base64 $i | sed ':a;N;$!ba;s/\n//g' >> $outDir/Data_and_Plots/Encoded/$my_file.txt
		echo -n "\" width=\"600\" height=\"600\" align=middle></iframe>" >> $outDir/Data_and_Plots/Encoded/$my_file.txt
	elif [[ $my_file =~ ".png" ]]; then # Encode PNGs
		echo -n "<img src=\"data:application/png;base64," > $outDir/Data_and_Plots/Encoded/$my_file.txt
		base64 $i | sed ':a;N;$!ba;s/\n//g' >> $outDir/Data_and_Plots/Encoded/$my_file.txt
		echo -n "\" />" >> $outDir/Data_and_Plots/Encoded/$my_file.txt
	fi
done

### Generate HTML file:
echo "
<!DOCTYPE html>
<html>
<body><h1>tsRNAsearch Report - Single Sample Analysis</h1>
<body><h2>${singleFile_basename}</h2>
<body><h3>tRNA Read Alignment Lengths</h3>
<br />
<embed src="$outDir/Data_and_Plots/${singleFile_basename}_tsRNA_tRNA-alignment-length.pdf" width="800px" height="800px" />
<br />
<body><h2>tRNA fragments (tsRNAs)</h2>
<br />
<body><h3>Distribution score</h3>
<br />
<embed src="$outDir/Data_and_Plots/${singleFile_basename}_tsRNA_Results.high-distribution-score.pdf" width="800px" height="800px" />
<br />
<body><h3>Cleavage score</h3>
<br />
<embed src="$outDir/Data_and_Plots/${singleFile_basename}_tsRNA_Results.high-cleavage-score.pdf" width="800px" height="800px" />
<br />
<body><h3>tsRNA Coverage Plots</h3>
<br />
<embed src="$outDir/Data_and_Plots/${singleFile_basename}_tsRNA_Coverage-plots.pdf" width="800px" height="800px" />
<br />
<body><h3>Predicted tsRNA type</h3>
$tsRNAtype
<br />
<body><h2>ncRNA fragments</h2>
<br />
<body><h3>Distribution score</h3>
<br />
<embed src="$outDir/Data_and_Plots/${singleFile_basename}_ncRNA_Results.high-distribution-score.pdf" width="800px" height="800px" />
<br />
<body><h3>Cleavage score</h3>
<br />
<embed src="$outDir/Data_and_Plots/${singleFile_basename}_ncRNA_Results.high-cleavage-score.pdf" width="800px" height="800px" />
<br />
<body><h3>ncRNA Coverage Plots</h3>
<br />
<br />Large file: $outDir/Data_and_Plots/${singleFile_basename}_ncRNA_Coverage-plots.pdf
<br />
</body>
</html>
" > $outDir/${singleFile_basename}.Results-summary.simple.html

echo -e "
---
title: 'tsRNAsearch Report'
output: 
  html_document:
    self_contained: no
    toc: true
    number_sections: true
    theme: cerulean
---

# Overview of tsRNAsearch analysis of ${singleFile_basename}

## tRNA Read Alignment Length

$(cat $outDir/Data_and_Plots/Encoded/${singleFile_basename}_tsRNA_tRNA-alignment-length.pdf.txt)

# tRNA fragments (tsRNAs)

## Distribution score

$(cat $outDir/Data_and_Plots/Encoded/${singleFile_basename}_tsRNA_Results.high-distribution-score.pdf.txt)

## Cleavage score

$(cat $outDir/Data_and_Plots/Encoded/${singleFile_basename}_tsRNA_Results.high-cleavage-score.pdf.txt)

## tsRNA Coverage Plots

$(cat $outDir/Data_and_Plots/Encoded/${singleFile_basename}_tsRNA_Coverage-plots.pdf.txt)

## Predicted tsRNA Type

\x60\x60\x60{r Summary report, echo=FALSE, error=TRUE}
library(knitr)
my.table <- read.csv(\"$outDir/tRNA-alignment/tsRNAs-classified-by-type_clean.txt\", sep = '\t')
kable(my.table, caption = 'Predicted tsRNA type')
\x60\x60\x60

# ncRNA fragments

## Distribution score

$(cat $outDir/Data_and_Plots/Encoded/${singleFile_basename}_ncRNA_Results.high-distribution-score.pdf.txt)

## Cleavage score

$(cat $outDir/Data_and_Plots/Encoded/${singleFile_basename}_ncRNA_Results.high-cleavage-score.pdf.txt)

## ncRNA Coverage Plots

All ncRNA coverage plots can be found here: 

* $outDir/Data_and_Plots/${singleFile_basename}_ncRNA_Coverage-plots.pdf

" > $outDir/Data_and_Plots/Run-Summary.Rmd
Rscript bin/Rmarkdown-to-HTML.R \
	$outDir/Data_and_Plots/Run-Summary.Rmd \
	$outDir/${singleFile_basename}.Results-summary.html

pdfunite \
	$outDir/Data_and_Plots/${singleFile_basename}_tsRNA_tRNA-alignment-length.pdf \
	$outDir/Data_and_Plots/${singleFile_basename}_tsRNA_Results.high-distribution-score.pdf \
	$outDir/Data_and_Plots/${singleFile_basename}_tsRNA_Results.high-cleavage-score.pdf \
	$outDir/Data_and_Plots/${singleFile_basename}_tsRNA_Coverage-plots.pdf \
	$outDir/Data_and_Plots/${singleFile_basename}_ncRNA_Results.high-distribution-score.pdf \
	$outDir/Data_and_Plots/${singleFile_basename}_ncRNA_Results.high-cleavage-score.pdf \
	$outDir/${singleFile_basename}.Results-summary.simple.pdf # This is the output file 

echo "Finished analysing "$singleFile" on $(date)" # Print pipeline end-time
echo "_____________________________________________________________________________________________________________________

"

