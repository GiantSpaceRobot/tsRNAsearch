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
	mkdir -p temp_files
	file_base=$(basename $f)
	filename=${file_base%.*}
	#echo $filename
	#./tiRNA-pipeline.sh -s $f -o Results/$filename -T
    #wait
	cat Results/$filename/HTSeq-count-output/*.count | grep -v ^__ | sort -k1,1 > temp_files/$filename.all_features.count
	sed -i '1s/^/Features\t'"$filename"'\n/' temp_files/$filename.all_features.count
	readsMapped=$(awk '{sum+=$2} END{print sum;}' temp_files/$filename.all_features.count)
	echo $readsMapped >> Results/$filename/Stats.log
done

### Gather count files
awk '{print $1}' temp_files/$filename.all_features.count > temp_files/HTSeq.all_features
for f in temp_files/*count; do
	head $f
	awk '{print $2}' $f | paste temp_files/HTSeq.all_features - >> temp_files/HTSeq.temp
	mv temp_files/HTSeq.temp temp_files/HTSeq.all_features
done

mv temp_files/HTSeq.all_features Results/HTSeq.all_features.count
#rm -rf temp_files

echo "Finished at $(date)"

