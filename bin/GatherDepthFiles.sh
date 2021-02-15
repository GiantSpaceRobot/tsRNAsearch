#!/bin/bash

count=0
for fname in $(cat $1); do
	count=$((count + 1))
	cond1="$(cut -d ',' -f1 <<< $fname)"
	substr=""
	case $cond1 in *.fastq.gz) cond1base=$(echo "${cond1/.fastq.gz/$substr}");; esac # If fastq.gz suffix in filename, remove
	case $cond1 in *.fq.gz) cond1base=$(echo "${cond1/.fq.gz/$substr}");; esac # If fq.gz suffix in filename, remove
	case $cond1 in *.fastq) cond1base=$(echo "${cond1/.fastq/$substr}");; esac # If fastq suffix in filename, remove
	case $cond1 in *.fq) cond1base=$(echo "${cond1/.fq/$substr}");; esac # If fq suffix in filename, remove
	case $cond1 in *.fa) cond1base=$(echo "${cond1/.fa/$substr}");; esac # If fq suffix in filename, remove
	case $cond1 in *.fa.gz) cond1base=$(echo "${cond1/.fa.gz/$substr}");; esac # If fq suffix in filename, remove
	echo "Finding file for $fname ($cond1base)..."
	find -L . -type f -name "${cond1base}*tRNAs_sorted.depth" -exec cp {} ${2}_file${count}_tsRNA_depth.readspermil \; & # Gather tsRNAs
	find -L . -type f -name "${cond1base}*ncRNAs_sorted.depth" -exec cp {} ${2}_file${count}_ncRNA_depth.readspermil \; & # Gather ncRNAs
	wait
	#mapped=$(grep "mapped" $outDir/$cond1base/Stats.log | awk '{print $3}')
	#mv ${2}_file${count}_tsRNA.depth ${2}_file${count}_tsRNA_depth.readspermil
	#mv ${2}_file${count}_ncRNA.depth ${2}_file${count}_ncRNA_depth.readspermil
	### Left-over reads:
	cp ${cond1base}_trimmed_tRNAs-almost-mapped_sorted.depth ${2}_${cond1base}_tRNAs-almost-mapped_RPM.depth
done
wait
### Concatenate data horizontally
paste ${2}_file*_tsRNA_depth.readspermil > tsRNA_${2}_concatenated.depth &
paste ${2}_file*_ncRNA_depth.readspermil > ncRNA_${2}_concatenated.depth &
paste ${2}_*_tRNAs-almost-mapped_RPM.depth > Combined_${2}_tRNAs-almost-mapped_RPM.depth & # leftover reads

cat ${2}_file*_tsRNA_depth.readspermil | sort -k1,1 -k2,2n > tsRNA_${2}_concatenated.depthVert &
cat ${2}_file*_ncRNA_depth.readspermil | sort -k1,1 -k2,2n > ncRNA_${2}_concatenated.depthVert &
wait
### Leftovers:
MeanCalculator.py \
	Combined_${2}_tRNAs-almost-mapped_RPM.depth \
	Combined_${2}_tRNAs-almost-mapped_RPM_depth.mean
### Calculate mean
MeanCalculator.py \
	tsRNA_${2}_concatenated.depth \
	tsRNA_${2}_concatenated_depth.mean &
MeanCalculator.py \
	ncRNA_${2}_concatenated.depth \
	ncRNA_${2}_concatenated_depth.mean &
wait
### Sort output
sort -k1,1 -k2,2n \
	tsRNA_${2}_concatenated_depth.mean > sorted_${2}_concatenated_tsRNA_depth.mean &
sort -k1,1 -k2,2n \
	ncRNA_${2}_concatenated_depth.mean > sorted_${2}_concatenated_ncRNA_depth.mean &
sort -k1,1 -k2,2n \
	Combined_${2}_tRNAs-almost-mapped_RPM_depth.mean > sorted_${2}_tRNAs-almost-mapped_depth.mean
wait

### Get mean and standard deviation
# tRNAs
Mean_Stdev.R \
	tsRNA_${2}_concatenated.depth \
	${2}_concatenated_mean_stdev_tsRNA.depth &
Mean_Stdev.R \
	Combined_${2}_tRNAs-almost-mapped_RPM.depth \
	Multimappers_${2}_concatenated_mean_stdev_tsRNA.depth &
# ncRNAs
Mean_Stdev.R \
	ncRNA_${2}_concatenated.depth \
	${2}_concatenated_mean_stdev_ncRNA.depth &
wait

### Create files with all RPM values for every coordinate of every feature in each condition
cat \
	Multimappers_${2}_concatenated_mean_stdev_tsRNA.depth \
	${2}_concatenated_mean_stdev_tsRNA.depth \
	| sort -k1,1 -k2,2n \
	> Everything_${2}_concatenated_mean_stdev_tsRNA.depth
cat \
	Everything_${2}_concatenated_mean_stdev_tsRNA.depth \
	${2}_concatenated_mean_stdev_ncRNA.depth \
	> Everything_${2}.depth
MeanCalculator.py \
	Everything_${2}.depth \
	Everything_${2}_depth.mean

### Sort output
sort -k1,1 -k2,2n \
	Everything_${2}_depth.mean > sorted_Everything_${2}_depth.mean &
sort -k1,1 -k2,2n \
	Everything_${2}_depth.mean > sorted_Everything_${2}_depth.mean &
wait

### Split sorted Everything files into tRNAs and other ncRNAs
# ncRNAs
grep ^ENS sorted_Everything_${2}_depth.mean \
	> sorted_Everything_ncRNAs_${2}_depth.mean
# tRNAs
grep -v ^ENS sorted_Everything_${2}_depth.mean \
	> sorted_Everything_tRNAs_${2}_depth.mean

echo "Finished transforming depth files"

