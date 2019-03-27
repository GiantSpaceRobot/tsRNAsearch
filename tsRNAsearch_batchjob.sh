#!/bin/sh

#SBATCH -t 00:20:00
#SBATCH -N 1
#SBATCH -c 20
#SBATCH -A cslif032c
#SBATCH -o batch-output.txt
#SBATCH -e batch-error.txt
#SBATCH --mail-user=pauldonovan@rcsi.com
#SBATCH --mail-type=BEGIN,END
#SBATCH --job-name=tsRNAsearch

tsRNAsearch -s fastq.fq.gz -o outputDir -p 20  
