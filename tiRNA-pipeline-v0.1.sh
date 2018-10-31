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

while getopts ":hs:1:2:o:" o; do
    case "${o}" in
		h)
			asciiArt
			usage
			info
			exit
			;;
		s)
			singleFile="$OPTARG"
			;;
		1)
			file1="$OPTARG"
			;;
		2)
			file2="$OPTARG"
			;;
		#f)
		#	#files=${OPTARG}
        #    set -f # disable glob
        #    IFS=',' # split on space characters
		#	array=($OPTARG)
		#	;;
		o)
			outDir="$OPTARG"
			;;
		*)
            echo "Error in input parameters!"
			usage
			exit 1
            ;;
    esac
done
#shift $((OPTIND-1))

if [ -z "$singleFile" ]; then # If the singleEnd variable is empty
	pairedEnd="True"
else
	pairedEnd="False"
fi

echo "Started at $(date)" # Print pipeline start-time

if [ ! -d $outDir ]; then # Check if the output directory exists. If not, create it
	mkdir $outDir
	mkdir $outDir/trim_galore_output
fi

# Run Trim_Galore
if [[ $pairedEnd = "True" ]]; then 
	trim_galore -o $outDir/trim_galore_output/ --paired $file1 $file2  
elif [[ $pairedEnd = "False" ]]; then
	trim_galore -o $outDir/trim_galore_output/ $singleFile
fi


