#!/usr/bin/env python
 
"""
Bedgraph thing
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
 
import sys
import numpy as np

bedgr = open(sys.argv[1], "r")
newfile = open(sys.argv[2], "w")
#newfile.write("Feature\tMean\tStandard deviation\tCoefficient of variation\n")
readlines = bedgr.readlines()
bedgr.close()

tRNAs = dict()

for line in readlines:
    strpline = line.strip().split("\t")
    if strpline[0] in tRNAs:
        tRNAs[strpline[0]].append(strpline[2])
    else:
        tRNAs[strpline[0]] = [strpline[2]] 

for k, v in tRNAs.iteritems():
    v2 = v[10:]
    v3 = v2[:len(v2) - 10] # Get values without the counts for the flanking Ns from the FASTA
    tRNALen = len(v) - 20     # Minus 20 because the FASTA was buffered/flanked with 10 Ns either side
    values = map(int, v3)
    meanV = np.mean(values)
    stdV = np.std(values)
    if float(meanV) == 0:
        coef_variation = "NA"
    else:
        coef_variation = float(stdV)/float(meanV)
    newfile.write(k + "\t" + str(meanV) + "\t" + str(stdV) + "\t" + str(coef_variation) + "\n")	
newfile.close()

