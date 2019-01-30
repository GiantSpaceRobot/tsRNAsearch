#!/usr/bin/env python
 
"""
Calculate the standard deviation of the mean of read counts
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"

import sys

data1 = open(sys.argv[1], "r")
newfile = open(sys.argv[2], "w")
readlines = data1.readlines()
data1.close()

### Calculate % difference between mean of two conditions
### Added 0.1 to all values to circumvent division by 0
for line in readlines:
    items = line.strip().split("\t")
    feature, mean1, mean2 = items[0], float(items[2]) + 0.1, float(items[5]) + 0.1
    percent_diff = (float(mean1)/float(mean2)) * 100
    newfile.write(feature + "\t" + str(mean1) + "\t" + str(mean2) + "\t" + str(percent_diff) + "\n")
newfile.close()



