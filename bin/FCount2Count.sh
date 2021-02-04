#!/bin/bash

grep -v featureCounts $1 \
	| grep -v ^Geneid \
	| awk -v OFS="\t" '{print $1, $7}' \
	> $1.count
