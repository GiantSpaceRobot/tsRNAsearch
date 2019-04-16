#!/usr/bin/env python
 
"""
Calculate the relative difference of the mean of read counts
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
__version__ = "Version 4"

import sys

data1 = open(sys.argv[1], "r")
newfile = open(sys.argv[2], "w")
readlines = data1.readlines()
data1.close()

### Calculate relative difference between mean of two conditions
for line in readlines:
    items = line.strip().split("\t")
    feature, mean1, mean2, std1, std2 = items[0], float(items[2]), float(items[6]), float(items[3]), float(items[7])
    diff = mean1 - mean2
    if mean1 == float(0) or mean2 == float(0):
        rel_diff = "NA"
    elif diff == float(0):
        rel_diff = 0
    else:
        rel_diff = (diff)/(min(mean1, mean2)) * 100
    newfile.write(feature + "\t" + str(mean1) + "\t" + str(mean2) + "\t" + str(rel_diff) + "\t" + str(std1) + "\t" + str(std2) + "\n")
newfile.close()



