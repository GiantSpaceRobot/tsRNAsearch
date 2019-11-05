#!/usr/bin/env python2
 
"""
Collapse similar tRNAs into a single feature for plotting 
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
__version__ = "Version 2"

import sys

bedgr = open(sys.argv[1], "r")
newfile = open(sys.argv[2], "w")
readlines = bedgr.readlines()
bedgr.close()

tRNAs = dict()
collapsed = dict()

for line in readlines:
    strpline = line.strip().split("\t")
    if strpline[0] in tRNAs:
        tRNAs[strpline[0]].append(strpline[2])
    else:
        tRNAs[strpline[0]] = [strpline[2]] 
    # Define tRNA groups:
    tRNA = strpline[0].split("-")[1]
    coord = strpline[1]
    reads = strpline[2].split() # Add number of reads to a list on it's own
    if tRNA not in collapsed:
        collapsed[tRNA] = {}
        collapsed[tRNA][coord] = reads
    else:
        ## If the dict exists in dict, update value. Else add new key
        if coord not in collapsed[tRNA]:
            collapsed[tRNA][coord] = reads
        else:
            val = (collapsed[tRNA][coord])
            newval = val + reads
            collapsed[tRNA][coord] = newval

for k, v in collapsed.iteritems():
    for key, values in v.iteritems():
        float_values = map(float, values)
        summed = sum(float_values)
        newfile.write("%s\t%s\t%s\n" % (k, key, summed))
newfile.close()

