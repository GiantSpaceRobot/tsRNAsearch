#!/usr/bin/env python2
 
"""
Calculate the mean and standard dev for 
every feature in a bedgraph/samtools depth file 
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
 
import sys
import numpy as np

bedgr = open(sys.argv[1], "r")
newfile = open(sys.argv[2], "w")
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
    v3 = v
    tRNALen = len(v)
    values = map(float, v3)
    meanV = np.mean(values)
    stdV = np.std(values)
    if float(meanV) == 0:
        coef_variation = "NA"
    else:
        coef_variation = float(stdV)/float(meanV)
    newfile.write(k + "\t" + str(meanV) + "\t" + str(stdV) + "\t" + str(coef_variation) + "\n")	
newfile.close()

