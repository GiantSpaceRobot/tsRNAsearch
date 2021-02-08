#!/bin/bash

### Create new name for output
newname=$(echo $1 | awk -F '_trimmed_accepted_hits_tRNAs' '{print $1}')
echo Filename: $newname

# Get nunmber of mapped reads from txt file
mapped=$(grep "$newname" $2 | awk '{print $2}')
echo Mapped reads: $mapped

### Arg 1 = Count file for all tRNAs
### Arg 2 = File with total read counts
### Arg 3 = Output file
CountCollapse.py \
	$1 \
	$3 \
	$mapped

