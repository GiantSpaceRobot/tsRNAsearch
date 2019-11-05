#!/usr/bin/env python2
 
"""
Take multimappers depth file and output list of 
expressed features to plot and list of all features
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
 
import sys

multimappers = open(sys.argv[1], "r")
newfile = open(sys.argv[2], "w")
allfile = open(sys.argv[3], "w")
readlines = multimappers.readlines()
multimappers.close()

tRNAs = dict()

for line in readlines:
    strpline = line.strip().split("\t")
    if strpline[0] in tRNAs:
        tRNAs[strpline[0]].append(strpline[2])
    else:
        tRNAs[strpline[0]] = [strpline[2]] 

allSet = set()
topSet = set()

for k,v in tRNAs.iteritems():
    if all(x == v[0] for x in v) == True:
        allSet.add(k)
    else:
        topSet.add(k)
        allSet.add(k)

for i in topSet:
    newfile.write(i + "\n")

for i in allSet:
    allfile.write(i + "\n")

newfile.close()
allfile.close()
