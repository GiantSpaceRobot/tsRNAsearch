#!/bin/bash

# Usage:
# Author: Paul Donovan 
# Email: pauldonovan@rcsi.com
# 19-Oct-2018

usage() { echo "Usage: $0 -p seq_1.fq seq_2.fq " 1>&2; exit 1; }
#Usage could also print out a manual on how the pipeline functions (link to manual file inside print statement)

while getopts ":hp:" o; do
    case "${o}" in
		h)
			echo "$usage"
			exit
			;;
		p)
            p=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${p}" ]; then   #if the string is empty, print usage
    usage
fi

#echo "Started at $(date)"
#StartTime="Pipeline initiated at $(date)"

#trim_galore -o trim_galore_output/ --paired Reads/SRR5788497_1_10000-reads.fq Reads/SRR5788497_2_10000-reads.fq
