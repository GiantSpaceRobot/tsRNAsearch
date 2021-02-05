#!/bin/bash

newname=$(echo $1 | awk -F '_mean-std_sorted.tsv' '{print $1}')

cp $1 ${newname}_depth.inf
echo -e "Feature\tMean\tStandard Deviation\tCoefficient of Variation" \
	> Header.txt
cat Header.txt \
	${newname}_depth.inf \
	> ${newname}_depth_stats.tsv
