#!/bin/sh
# Usage: ./Master.sh /path/to/directory/with/files/
# Author: Paul Donovan 
# Email: pauldonovan@rcsi.com
# 7-12-2018

echo "Started at $(date)"
StartTime="Pipeline initiated at $(date)" 

### Precautionary measures before running analysis
if [ $# -eq 0 ]; then
	echo "No arguments provided"
	exit 1
fi

for f in $1/*; do
	file_base=$(basename $f)
	filename=${file_base%.*}
	echo $filename
	./tiRNA-pipeline.sh -s $f -o Results/$filename -T
    wait
done

echo "Finished at $(date)"

