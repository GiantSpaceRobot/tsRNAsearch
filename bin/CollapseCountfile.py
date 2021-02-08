#!/usr/bin/env python2
 
"""
Collapse tRNA counts for individual tRNA genes into isodecoders 
(e.g. collapse chr11.trna418-ArgCCG and chr13.trna488-ArgCCG into ArgCCG)
Output is for use with DESeq2
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"
 
import sys

countfile = open(sys.argv[1], "r")
countfile_readlines = countfile.readlines()
countfile.close()
newfile = open(sys.argv[2], "w")

collapse = dict()

for line in countfile_readlines:
    if line.startswith("Features"):
        pass
    else:
        strpline = line.strip().split("\t")
        read_count = float(strpline[1])
        if line.startswith("chr"):
            isodecoder = strpline[0].split("-")[1]
            if isodecoder in collapse:
                previousCount = collapse.get(isodecoder)
                new_read_count = float(previousCount) + float(read_count) #Add all counts for this isodecoder together
                collapse[isodecoder] = new_read_count
            else:
                collapse[isodecoder] = read_count
        else:
            collapse[strpline[0]] = read_count
        
for key in sorted(collapse): #Sort and print the output
    newfile.write(key + "\t" + str(collapse[key]) + "\n")

newfile.close()

