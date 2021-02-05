#!/bin/bash

### Create base name for input
newname=$(echo $1 | awk -F '_Collapsed.bam' '{print $1}')

### Convert BAM to SAM
samtools view -h -o ${newname}.sam $1

echo Sample name: ${newname}

### Split SAM by gene name (genes beginning in 'ENS' are ncRNAs, the rest are tRNAs)
echo "Split SAM by gene name (genes beginning in 'ENS' are ncRNAs, the rest are tRNAs)"
grep ^@ ${newname}.sam | grep -v 'chr.*trna' > SamHeader.sam &
grep ENS ${newname}.sam | grep -v 'chr.*trna' | grep -v ^@ > ncRNAs.sam &
grep 'chr.*trna' ${newname}.sam > tsRNAs_aligned.sam &
wait
cat SamHeader.sam ncRNAs.sam \
	> ncRNAs_aligned.sam
### Generate BAM from tRNA SAM
echo "Generate BAM from tRNA SAM"
samtools view -bS tsRNAs_aligned.sam \
	| samtools sort \
	> ${newname}_accepted_hits_tRNAs.bam
#samtools index accepted_hits_tRNAs.bam &
### Generate BAM from ncRNA SAM
echo "Generate BAM from ncRNA SAM"
samtools view -bS ncRNAs_aligned.sam \
	| samtools sort \
	> ${newname}_accepted_hits_ncRNAs.bam
#samtools index accepted_hits_ncRNAs.bam &
echo "Remove intermediate files"
rm \
	SamHeader.sam \
	ncRNAs.sam \
	tsRNAs_aligned.sam \
	ncRNAs_aligned.sam &
