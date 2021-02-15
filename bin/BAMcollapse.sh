#!/bin/bash

### Create new name for output
newname=$(echo $1 | cut -f 1 -d '.')

### Convert BAM to SAM
samtools view -h -o ${newname}.sam $1

### How many chunks to create:
readNumber=$(grep -v ^@ ${newname}.sam | awk '{print $1}' | uniq | wc -l)
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
cores=$2
if (( $cores > $chunks )); then
	# make sure cores is not set higher than no. of files after split
	cores=$chunks
fi
echo "Using $cores cores"
fileLen=$(< "${newname}.sam" wc -l)
division1=$((fileLen/chunks))
division=$((division1 + 1))
echo "Splitting SAM file into $division chunks"
myFile="tempFile"
mkdir -p tempDir

### Remove header
grep ^@ ${newname}.sam > tempDir/myHeader.txt &
grep -v ^@ ${newname}.sam > tempDir/mySAM.sam &
wait

### 
fileToCollapse=tempDir/mySAM.sam

### Split file
split -l $division $fileToCollapse tempDir/splitFile_

### Gather first and last read from every split file and add to separate file. Remove these reads from the split files.
echo "Gather first and last read from every split file and add to separate file."
for i in tempDir/splitFile_*; do
	base=$(basename $i)
	first=$(awk 'NR==1' $i | awk '{print $1}') 
	echo $first >> tempDir/${myFile}_HeadsAndTails.txt
	last=$(awk 'END{print}' $i | awk '{print $1}') 
	echo $last >> tempDir/${myFile}_HeadsAndTails.txt
done
sort tempDir/${myFile}_HeadsAndTails.txt | uniq > tempDir/${myFile}_HeadsAndTails_uniq.txt #remove duplicates
sed -i 's/$/\t/' tempDir/${myFile}_HeadsAndTails_uniq.txt # Add tab to end of every line to match pattern exactly
grep -f tempDir/${myFile}_HeadsAndTails_uniq.txt $fileToCollapse > tempDir/edit_heads-and-tails #grep all patterns from the heads/tails file
for i in tempDir/splitFile_*; do
	base=$(basename $i)
	grep -v -f tempDir/${myFile}_HeadsAndTails_uniq.txt $i > tempDir/edit_${base}
done
linesInHeadTailFile=$(wc -l tempDir/edit_${base})
echo "There are $linesInHeadTailFile lines in the head & tail file
"
### Run SAMcollapse2.py. This loop will only run $cores processes at once
rm -f collapsed-reads.txt # Overwrite this file if it exists
COUNTER=0
for i in tempDir/edit_*; 
do
	base=$(basename $i)
	SAMcollapse.py $i ${fileToCollapse}_${base} >> collapsed-reads.txt  & 
	numjobs=($(jobs | wc -l))
	echo Running job number ${COUNTER} of ${chunks}... 
	COUNTER=$[$COUNTER + 1]
	#echo Initialising job number $numjobs
	#echo There are ${#background[@]} jobs now
    while (( $numjobs > $cores )); do
    	echo There are $numjobs jobs now. Waiting for jobs to finish...
		numjobs=($(jobs | wc -l))
		sleep 2 #Enter next loop iteration
	done
done
wait
readsCollapsedSpecies=$(awk '{split($0,a," "); sum += a[1]} END {print sum}' collapsed-reads.txt)
readsCollapsedGroup=$(awk '{split($0,a," "); sum += a[2]} END {print sum}' collapsed-reads.txt)
echo -e "SAM collapse results:\n\t$readsCollapsedSpecies reads collapsed at the tRNA species level (e.g. 2 gene copies of ProCCG)\n\t$readsCollapsedGroup reads collapsed at the tRNA group level (e.g. ProCCG and ProAAG)"
### Concatenate results
echo "Gathering reads that were mapped to similar tRNAs..."
echo -e "tRNA.group\tread.start\tread.end.approx\tread.name" \
	> ${newname}_tRNAs-almost-mapped.txt
cat tempDir/*tRNAs-almost-mapped* | sort \
	>> ${newname}_tRNAs-almost-mapped.txt
grep -v "tRNA.group" ${newname}_tRNAs-almost-mapped.txt \
	| awk '{print $1}' \
	| uniq -c \
	| awk '{print $2"\t"$1}' \
	> ${newname}_tRNAs-almost-mapped.count

mkdir tempDir/tRNAsAlmostMapped
mv tempDir/*tRNAs-almost-mapped* tempDir/tRNAsAlmostMapped/
echo "Concatenating SAM header with collapsed files..."
cat tempDir/myHeader.txt ${fileToCollapse}*_edit_* | samtools sort > ${newname}_Collapsed.bam
### Sort SAM to BAM
#samtools view -bS *Collapsed.sam > "$reads.simpleName"_Collapsed.bam
#samtools sort Collapsed.sam > Collapsed.bam

rm -rf tempDir/
echo "Finished collapsing SAM file"

