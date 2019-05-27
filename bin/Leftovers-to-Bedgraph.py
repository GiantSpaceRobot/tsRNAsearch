#!/usr/bin/env python
 
"""
Expand leftovers tRNA group files into Bedgraph format
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
__version__ = "Version 1"

import sys

leftovers = open(sys.argv[1], "r")
tRNAlengths = open(sys.argv[2], "r")
newfile = open(sys.argv[3], "w")
readlines = leftovers.readlines()
leftovers.close()

tRNAs = dict()
tRNALen = dict()
negSet = set()

### Create dict from known tRNA length file
for i in tRNAlengths:
    ispl = i.strip().split("\t")
    tRNAgroup, tRNAlen = ispl[0], ispl[1]
    tRNALen[tRNAgroup] = tRNAlen
    negSet.add(tRNAgroup)

### Create read counts to dict for every tRNA group (e.g. Val: [1,1,2,2,3,3,4,5])
for line in readlines:
    strpline = line.strip().split("\t")
    if strpline[0].startswith("tRNA.group"):
        pass
    else:
        tRNA = strpline[0]
        start = strpline[1]
        stop = strpline[2]
        read = strpline[3]
        read_coverage = range(int(start), int(stop) + 1)
        read_cov_list = list()
        read_cov_list.append(read_coverage)
        if tRNA in tRNAs:
            values = tRNAs.get(tRNA)
            values.extend(read_coverage)
            tRNAs[tRNA] = values
        else:
            tRNAs[tRNA] = read_coverage 

### Write depth results for all tRNAs with reads multi-mapped
for k,v in tRNAs.iteritems():
    negSet.remove(k) #remove tRNAs that have reads mapped from the negSet
    getLen = int(tRNALen.get(k))  # Get the total length of each tRNA using a prebuilt tRNA length file
    myRange = range(1,getLen+1)
    for i in myRange:
        count = v.count(i)
        newfile.write("%s\t%s\t%s\n" % (k, i, count))
### Write empty depth results for all tRNAs that had no multi-mapping reads
### This step is important because it ensures all resulting depth files are the same length
for negtRNA in negSet:
    getLen = int(tRNALen.get(negtRNA))  # Get the total length of each tRNA using a prebuilt tRNA length file
    myRange = range(1,getLen+1)
    for i in myRange:
        newfile.write("%s\t%s\t%s\n" % (negtRNA, i, 0))
newfile.close()

