#!/usr/bin/env python
 
"""
Calculate the mean of read counts for concatenated genomecov files
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"

import sys

data1 = open(sys.argv[1], "r")
newfile = open(sys.argv[2], "w")
readlines = data1.readlines()
readlines2 = readlines
data1.close()

feature_set = set(line.strip().split("\t")[0] for line in open(sys.argv[1]))

for feature in feature_set:
    feature_length = 0
    for line in readlines:
        if feature in line:
            feature_length = feature_length + 1
    for line in readlines:
        if feature in line:
            items = line.strip().split("\t")
            if int(items[1]) > 10:
                if int(items[1]) <= int(feature_length - 10):
                    #print ("Feature: %s\nLine of first occurence: %s\nPrinted line: %sFeature length: %s\n\n" % (feature, linecount, firstline, feature_length))
                    #items = line.strip().split("\t")
                    data = items[2::3]
                    data = map(float, data)
                    summed = sum(data)
                    mean = (float(summed) / len(data))
                    newfile.write(items[0] + "\t" + str(int(items[1]) - 10) + "\t" + str(mean) + "\n")

newfile.close()
