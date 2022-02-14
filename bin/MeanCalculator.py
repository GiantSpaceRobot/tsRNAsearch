#!/usr/bin/env python2
 
"""
Calculate the mean of read counts for concatenated samtools depth files
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
__version__ = "Version 5"

import sys, numpy

data1 = open(sys.argv[1], "r")
newfile = open(sys.argv[2], "w")
readlines = data1.readlines()
readlines2 = readlines
data1.close()

for line in readlines:
    items = line.strip().split("\t")
    feature = items[0]
    data = items[2::3]
    data = [0 if not x.isdigit() else x for x in data] # convert NA or other nan to zero
    data = map(float, data)
    summed = sum(data)
    mean = (float(summed) / len(data))
    ### Standard deviation:
    stdev = numpy.std(data)
    newfile.write(items[0] + "\t" + str(items[1]) + "\t" + str(mean) + "\t" + str(stdev) + "\n")

newfile.close()
