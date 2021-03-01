#!/bin/bash

### Normalise depth file by transcripts per kilobase million (TPM)
#arg1 depthfile
#arg2 simple name
#arg3 tRNA gtf
#arg4 ncRNA gtf

# Create new name for output
newname=$(echo $1 | awk -F '_trimmed' '{print $1}')
echo Filename: $newname

# Get nunmber of mapped reads from txt file
#mapped=$(grep "$newname" $2 | awk '{print $2}')
#echo Mapped reads: $mapped

# Convert depth to RPM
Depth-to-TPM.R \
	$1 \
	$2 \
	$3 \
	$4

echo "Finished normalising depth file"


