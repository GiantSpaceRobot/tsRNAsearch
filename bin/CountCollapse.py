#!/usr/bin/env python2
 
"""
Collapse count file to amino acid species
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
__version__ = "Version 1"

import sys, numpy

data1 = open(sys.argv[1], "r")
newfile = open(sys.argv[2], "w")
readlines = data1.readlines()
data1.close()

my_dict = dict()
tRNA_mapped_reads = 0
total_mapped_reads = int(sys.argv[3])

for line in readlines:
    items = line.strip().split("\t")
    feature = items[0]
    tRNA_mapped_reads = tRNA_mapped_reads + int(items[1])
    if feature.startswith("chr"): # If feature is a tRNA
        if "iMet" in feature: # If feature is iMethionine
            feature = feature.split("-")[1][:4]
        else:
            feature = feature.split("-")[1][:3]
    if feature in my_dict:
        new_total = my_dict[feature] + int(items[1])
        my_dict[feature] = new_total
    else:
        my_dict[feature] = int(items[1])

# Create output file
newfile.write("tRNA\traw.read.count\tpercent.as.tRNA.total\tpercent.as.mapped.read.total\n")
readcount = 0
total_tRNApercent = 0
total_totalpercent = 0
for k,v in my_dict.iteritems():
    percent_of_tRNAs = round((float(v)/float(tRNA_mapped_reads)) * 100, 2)
    percent_of_total = round((float(v)/float(total_mapped_reads)) * 100, 2)
    newfile.write(k + "\t" + str(v) + "\t" + str(percent_of_tRNAs) + "\t" + str(percent_of_total) + "\n")    
    readcount = readcount + int(v)
    total_tRNApercent = total_tRNApercent + float(percent_of_tRNAs)
    total_totalpercent = total_totalpercent + float(percent_of_total)
newfile.write("Total\t" + str(readcount) + "\t" + str(total_tRNApercent) + "\t" + str(total_totalpercent) + "\n")
newfile.close()
