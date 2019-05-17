# Collapses SAM files

import sys

SAMList = list(line.strip().split("\t") for line in open(sys.argv[1]))
Output = open(sys.argv[2], "w")

myDict = dict()

for samLine in SAMList:
    line = "\t".join(samLine)
    if samLine[0].startswith("@"): #if a header
        Output.write(line + "\n")
    else:
        read = str(samLine[0])#.split("_simread")[0])
        if read in myDict.keys():
            myDict[read].append(line)
        else:
            myList = [line]
            myDict[read] = (myList) 

tRNAgroupDict = dict()
counter = 0
counterGroup = 0
for k,v in myDict.iteritems():
    if v[0].split("\t")[2].startswith("chr"):
        feature = v[0].split("\t")[2].split("-")[1]
    else:
        feature = v[0].split("\t")[2]
    database_matches = list()
    for i in v:
        match = i.split("\t")[2]
        database_matches.append(match)
    if all(feature in x for x in database_matches) == True:  # if every value in the dict contains the feature species
        newAlignment = v[0].split(":")
        newAlignment[-1:] = ["1"] # Replace last element (NH flag DB match #) with 1
        newAlignment = ":".join(newAlignment)
        Output.write(newAlignment + "\n")     
        counter = counter + 1
    else:
        if v[0].split("\t")[2].startswith("chr"):
            tRNAgroup = (v[0].split("\t")[2].split("-")[1])[:3] # e.g. 'Pro' from a Proline tRNA
        else:
            tRNAgroup = v[0].split("\t")[2]
        if all(tRNAgroup in x for x in database_matches) == True:  # if every value in the dict contains the tRNA group
            counterGroup = counterGroup + 1
            newLine = v[0].split("\t")
            newLine = tRNAgroup + "\t" + newLine[3] + "\t" + str(len(newLine[9]) + int(newLine[3])) + "\t" + newLine[0]
            if tRNAgroup in tRNAgroupDict.keys():
                tRNAgroupDict[tRNAgroup].append(newLine)
            else:
                myList = [newLine]
                tRNAgroupDict[tRNAgroup] = (myList)
        #for line in v:
        #    Output.write(line + "\n")

Output.close()

print "SAMcollapse.py results:\n%s reads collapsed at the tRNA species level (e.g. 2 gene copies of ProCCG)\n%s reads collapsed at the tRNA group level (e.g. ProCCG and ProAAG)" % (counter, counterGroup)

leftovers = open(sys.argv[2] + "_tRNAs-almost-mapped.txt", "w")
#leftovers.write("tRNA.group\tread.start\tread.end.approx\tread.name\n")
for k,v in tRNAgroupDict.iteritems():
    for values in v:
        leftovers.write(values + "\n")
leftovers.close()

