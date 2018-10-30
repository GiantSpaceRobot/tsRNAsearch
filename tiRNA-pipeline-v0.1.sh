#!/bin/bash

# Usage:
# Author: Paul Donovan 
# Email: pauldonovan@rcsi.com
# 19-Oct-2018

asciiArt() { echo "
  __  .____________  _______      _____    __________.__              .__  .__               
_/  |_|__\______   \ \      \    /  _  \   \______   \__|_____   ____ |  | |__| ____   ____  
\   __\  ||       _/ /   |   \  /  /_\  \   |     ___/  \____ \_/ __ \|  | |  |/    \_/ __ \ 
 |  | |  ||    |   \/    |    \/    |    \  |    |   |  |  |_> >  ___/|  |_|  |   |  \  ___/ 
 |__| |__||____|_  /\____|__  /\____|__  /  |____|   |__|   __/ \___  >____/__|___|  /\___  >
                 \/         \/         \/               |__|        \/             \/     \/ 
" 1>&1; }
usage() { echo "Usage: $0 -f seq_1.fq,seq_2.fq -o OutputDirectory" 1>&2; }
info() { echo "
Options
	
	-h	Print the usage and options information
	-f	File(s) for analysis. Please use comma-delimited format for paired-end files
	-o	Output directory for the results and log files
	" 1>&2; }

while getopts ":hf:o:" o; do
    case "${o}" in
		h)
			asciiArt
			usage
			info
			exit
			;;
		f)
			files=${OPTARG}
            set -f # disable glob
            IFS=',' # split on space characters
            array=($OPTARG)
			;;
		o)
			outDir={OPTARG}
			;;
		*)
            echo "Error in input parameters!"
			usage
			exit 1
            ;;
    esac
done
shift $((OPTIND-1))

pairedEnd = True

if [ ${#array[@]} -gt 2 ]; then
	echo "Error, ${#array[@]} files supplied. Only 1 (single-end) or 2 (paired-end data) files accepted."
	usage
	exit
elif [ ${#array[@]} -lt 1 ]; then  
    echo "Error, no files supplied"
	usage
	exit
elif [ ${#array[@]} -eq 1 ]; then
	pairedEnd = False
fi

echo "Started at $(date)"
#StartTime="Pipeline initiated at $(date)"

#trim_galore -o trim_galore_output/ --paired Reads/SRR5788497_1_10000-reads.fq Reads/SRR5788497_2_10000-reads.fq
