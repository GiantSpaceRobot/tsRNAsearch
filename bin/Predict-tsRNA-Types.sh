#!/bin/bash

for i in $(cat $1)
do
	echo $i
	my_sample=$(echo $i | awk -F '_accepted_hits' '{print $1}')
	echo $my_sample
	grep -v ^feature ${my_sample}_tRNAs-almost-mapped_sorted.depth \
		> ${my_sample}_multimappers_noheader.depth
	cat $i \
		${my_sample}_multimappers_noheader.depth \
		> ${my_sample}_all-tRNAs.depth
	tsRNA-type-classification.R \
		${my_sample}_all-tRNAs.depth \
		${my_sample}_tsRNAs-classified-by-type.tsv
	awk '{print $1"\t"$2}' ${my_sample}_tsRNAs-classified-by-type.tsv \
		> ${my_sample}_tsRNAs-classified-by-type_clean.tsv
	sed -e 's/^/<br \/>/' ${my_sample}_tsRNAs-classified-by-type_clean.tsv \
		> ${my_sample}_tsRNAs-classified-by-type_clean_HTML.txt
	sed -e 's/^/* /' ${my_sample}_tsRNAs-classified-by-type_clean.tsv \
		> ${my_sample}_tsRNAs-classified-by-type_clean_RmdHTML.txt
done

