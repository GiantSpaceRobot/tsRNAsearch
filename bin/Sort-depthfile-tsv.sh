#!/bin/bash

### Create base name for input
newname=$(echo $1 | awk -F '.tsv' '{print $1}')

awk '$2>0' $1 \
    > ${newname}_mean-std.tsv
sort -k4,4nr ${newname}_mean-std.tsv \
    > $2

