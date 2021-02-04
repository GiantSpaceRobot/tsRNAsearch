#!/usr/bin/env python3
 
"""
Count tsRNAs from SAM
"""
 
__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovandonegal@gmail.com"
 
import sys
import re
import collections

# Function that outputs dictionary of tRNA_start-stop (e.g. Val_1-20)
# where the start stop must be 16 nt or longer
def tRNA_start_stop(tRNA, start, stop):
    """Given tRNA name and first and last nucleotide index, 
    output all possible tsRNAs 16 nt or greater"""
    start = int(start)
    stop = int(stop)
    tsRNA_dict = dict()
    for kmer_len in range(16, stop + 1):  # Iterate over kmer length from 16 to tRNA length (i.e. stop variable)
        for kmer_start in range(1, stop + 10 - kmer_len): # Go from first to last nucleotide (minus kmer length)
            ### I added stop + 10 because the FASTA and GTF were in disagreement. 
            ### This is not a problem as fake dict accounts will remain at value of 0
            possible_tsRNA = tRNA + "_" + str(kmer_start) + "-" + str(kmer_start + kmer_len)
            if possible_tsRNA not in tsRNA_dict:
                tsRNA_dict[possible_tsRNA] = 0
    return(tsRNA_dict)

giant_tsRNA_dict = dict()

# GTF file to create dictionary of possible tsRNAs
with open(sys.argv[1]) as f: 
    for line in f:
        spltline = line.strip().split("\t")
        current_tRNA = spltline[0].split("-")[1]
        current_tRNA_generic = current_tRNA[:3] # Get Gly if tRNA is GlyGCC
        tRNA_start_stop_output = tRNA_start_stop(current_tRNA, spltline[3], spltline[4]) # Run function on current tRNA from GTF
        giant_tsRNA_dict.update(tRNA_start_stop_output) # add new tsRNA_start_stops to big dict. Overwrite keys by default
        tRNAgeneric_start_stop_output = tRNA_start_stop(current_tRNA_generic, spltline[3], spltline[4]) # Run function on current tRNA from GTF
        giant_tsRNA_dict.update(tRNAgeneric_start_stop_output) # add new tsRNA_start_stops to big dict. Overwrite keys by default

for line in sys.stdin: 
    #for line in f:
    if line.startswith("@"): # This is a header. Skip
        pass
    else:
        splitline = line.strip().split("\t")
        tRNA = splitline[2].split("-")[1]
        mapping_from = int(splitline[3])
        cigar = splitline[5]
        #cigar = "12M10S18M4I5D17M" # Test string
        split_cigar = re.split('([0-9]+[S M D I H N]+)', cigar) # Split CIGAR string by AnyNumberAnyLetter pattern
        split_cigar = [i for i in split_cigar if i] # If element in list is real (i.e. not whitespace), keep it
        
        ### Remove S from start and end. They are read overhangs
        if "S" in split_cigar[0]: # If the first element contains S
            del split_cigar[0] # Delete this element
        if "S" in split_cigar[-1]:
            del split_cigar[-1]

        # Sum MSIN CIGARs and ignore I
        read_span_of_tRNA = 0
        for i in split_cigar:
            section_len = i[:-1]
            if "M" in i or "S" in i or "D" in i or "N" in i:
                read_span_of_tRNA = read_span_of_tRNA + int(section_len)
            elif "I" in i: # Do nothing
                pass
            else:
                print("ERROR: unknown CIGAR character: ", i)
                print("SAM account:\n", splitline)
        mapping_to = mapping_from + int(read_span_of_tRNA) - 1 # -1 to correct for mapping start position including the nucleotide it starts with (i.e. starting from 1 INCLUDES 1 in the read mapping)
        SAM_tsRNA = tRNA + "_" + str(mapping_from) + "-" + str(mapping_to)
        if SAM_tsRNA not in giant_tsRNA_dict:
            print("ERROR: ", SAM_tsRNA, " not in dictionary!")
            print("SAM account:\n", splitline)
        else:
            old_val = giant_tsRNA_dict.get(SAM_tsRNA)
            giant_tsRNA_dict[SAM_tsRNA] = int(old_val) + 1  # Add 1 to the dictionary counter

### Sort dict, write to file
sorted_dict = collections.OrderedDict(sorted(giant_tsRNA_dict.items())) # Sort dict by key
output = open(sys.argv[2], 'w')
for k,v in sorted_dict.items():
    output.write(k + "\t" + str(v) + "\n")
output.close()
