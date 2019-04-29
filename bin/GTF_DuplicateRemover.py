#!/usr/bin/env python
 
"""
Remove features from GTF that are found in provided text file
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
 
import sys

d = {}
with open(sys.argv[1]) as f:
    for line in f:
       key = line.strip()
       d[key] = 0

gff = open(sys.argv[2], "r")
newgff = open(sys.argv[3], "w")
readlines = gff.readlines()
gff.close()

for line in readlines:
    if line.startswith("#"):
        newgff.write(line)
    else:
        splitline = line.strip().split("\t")
        featureID = splitline[8].split('"')[1]
        if featureID in d:
            pass
            #print (featureID)
        else:
            #gff = (featureID[1], splitline[1], splitline[2], "1", str((int(splitline[4]) + 1) - int(splitline[3])), splitline[5], splitline[6], splitline[7], splitline[8])
            #newgff.write("\t".join(gff) + "\n")
            newgff.write(line)
newgff.close()
