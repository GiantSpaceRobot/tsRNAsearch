#!/bin/sh
# Usage: ./FindFungi-0.23.sh /path/to/FASTA-file.fastq File-name
# Author: Paul Donovan 
# Email: pauldonovandonegal@gmail.com
# 11-Jul-2018

echo "Started at $(date)"
StartTime="Pipeline initiated at $(date)" 

### Precautionary measures before running analysis
if [ $# -eq 0 ]; then
	echo "No arguments provided"
	exit 1
fi

for f in $1/*; do
    singleFile_base=${$1##*/}    # Remove pathname
	singleFile_basename="$( cut -d '.' -f 1 <<< '$singleFile_base' )" # Get filename before first occurence of .	
	File=$(basename $f)
    ./tiRNA-pipeline.sh -s $1 -o Results/$File -T
    wait
done

echo "Finished at $(date)"

