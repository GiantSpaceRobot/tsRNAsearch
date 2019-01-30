#!/usr/bin/env python
 
"""
Normalise Genomecov to RPM
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
 
import sys

htseq = open(sys.argv[1], "r")
htseq_readlines = htseq.readlines()
htseq.close()

#reads_mapped = open(sys.argv[2], "r")
#reads = reads_mapped.readlines()
#total_line = reads[3] # Get line with total number of reads mapped
#total_reads = float(total_line.strip().split("\t")[1])
total_reads = sys.argv[2]
total_reads = float(total_reads.strip().split("\n")[0])
scaling_factor = float(total_reads/1000000)
#print (scaling_factor)
#reads_mapped.close()


newfile = open(sys.argv[3], "w")

for line in htseq_readlines:
    strpline = line.strip().split("\t")
    read_count = float(strpline[2])
    rpm = read_count/scaling_factor
    newfile.write(strpline[0] + "\t" + strpline[1] + "\t" + str(rpm) + "\n")	
newfile.close()

