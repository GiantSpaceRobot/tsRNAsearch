#!/bin/bash

### Normalise depth file by reads per million (RPM)

# Create new name for output
newname=$(echo $1 | awk -F '_trimmed_accepted_hits' '{print $1}')
echo Filename: $newname

# Get nunmber of mapped reads from txt file
mapped=$(grep "$newname" $2 | awk '{print $2}')
echo Mapped reads: $mapped

# Convert depth to RPM
Depth-to-Depth_RPM.py \
	$1 \
	$mapped \
	$3

echo "Finished normalising depth file"


