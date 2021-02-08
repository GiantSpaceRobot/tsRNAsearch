#!/usr/bin/env python2

# Collapses SAM files

import sys

partial_count_file = list(line.strip().split("\t") for line in open(sys.argv[1]))
empty_count_file = list(line.strip().split("\t") for line in open(sys.argv[2]))
filled_count_file = open(sys.argv[3], "w")

trna_set = set()

for line in empty_count_file:
    trna = line[0].split('-')[1] # Isolate tRNA species
    # Extract first 3 letters of tRNA name
    if trna.startswith("i"):
        trna = trna[:4]
    else:
        # If tRNA is iMet..., extract iMet
        trna = trna [:3]
    trna_set.add(trna)

trna_species_dict = dict()

# Add multimapper read counts to dict
for v in partial_count_file:
    key = v[0]
    value = v[1]
    trna_species_dict[key] = (value)    

# Fill the remainder of new dict with empty values
for i in trna_set:
    if i in trna_species_dict.keys():
        pass
    else:
        trna_species_dict[i] = (0)

# Write dict to file
for k,v in trna_species_dict.iteritems():
    filled_count_file.write(k + "\t" + str(v) + "\n")
filled_count_file.close()
