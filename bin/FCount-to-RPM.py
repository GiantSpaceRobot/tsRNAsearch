#!/usr/bin/env python
 
"""
Convert FeatureCounts into RPM
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
 
import sys

htseq = open(sys.argv[1], "r")
htseq_readlines = htseq.readlines()
htseq.close()

total_reads = float(sys.argv[2])
scaling_factor = float(total_reads/1000000)

raw_and_rpm_outfile = open(sys.argv[3] + ".raw_rpm.txt", "w")
rpm_outfile = open(sys.argv[3] + ".rpm.count", "w")

for line in htseq_readlines:
    strpline = line.strip().split("\t")
    if line.startswith("Feature"):
        raw_and_rpm_outfile.write(strpline[0] + "\t" + strpline[1] + "\tReadsPerMillion\n")
    else:
        read_count = int(strpline[1])
        if line.startswith("__"):
            pass
        else:
            rpm = read_count/scaling_factor
            raw_and_rpm_outfile.write(strpline[0] + "\t" + str(read_count) + "\t" + str(rpm) + "\n")	
            rpm_outfile.write(strpline[0] + "\t" + str(rpm) + "\n")	
raw_and_rpm_outfile.close()
rpm_outfile.close()
