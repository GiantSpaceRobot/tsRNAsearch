#!/bin/bash

#cat tRNA_Counts.txt | while read line;

for i in $(cat $1)
do
	echo $i
	my_sample=$(echo $i | awk -F '_accepted_hits' '{print $1}')
	echo $my_sample
	cat ${i} ${my_sample}_tRNAs-almost-mapped.count ${my_sample}_accepted_hits_ncRNAs.bam.fcount.count > ${my_sample}_all-counts.count
	reads_mapped=$(awk '{sum+=$2} END{print sum;}' ${my_sample}_all-counts.count)
	echo -e "${my_sample}\t${reads_mapped}" > ${my_sample}_reads-mapped.txt
done

