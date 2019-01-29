#!/usr/bin/env python
 
"""
Calculate the mean of read counts for concatenated genomecov files
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
__version__ = "Version 3"

import sys

data1 = open(sys.argv[1], "r")
newfile = open(sys.argv[2], "w")
readlines = data1.readlines()
readlines2 = readlines
data1.close()

d = {}
with open(sys.argv[1]) as f:
    for line in f:
        (key, val) = (line.split("\t")[0], line.split("\t")[1])
        if key in d:
            # append the new number to the existing key at this slot
            d[key].append(val)
        else:
            # create a new array in this slot
            d[key] = [val]

for line in readlines:
    items = line.strip().split("\t")
    feature = items[0]
    feature_length = (len(d.get(feature)))
    if int(items[1]) > 10:
        if int(items[1]) <= int(feature_length - 10):
            #print ("Feature: %s\nLine of first occurence: %s\nPrinted line: %sFeature length: %s\n\n" % (feature, linecount, firstline, feature_length))
            data = items[2::3]
            data = map(float, data)
            summed = sum(data)
            mean = (float(summed) / len(data))
            newfile.write(items[0] + "\t" + str(int(items[1]) - 10) + "\t" + str(mean) + "\n")

newfile.close()
