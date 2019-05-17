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

for i in tRNAlengths:
    ispl = i.strip().split("\t")
    tRNAgroup, tRNAlen = ispl[0], ispl[1]
    tRNALen[tRNAgroup] = tRNAlen

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

for k,v in tRNAs.iteritems():
    getLen = int(tRNALen.get(k))  # Get the total length of each tRNA using a prebuilt tRNA length file
    myRange = range(1,getLen+1)
    for i in myRange:
        count = v.count(i)
        newfile.write("%s\t%s\t%s\n" % (k, i, count))
        #print (k, i, count)
    #break
newfile.close()

