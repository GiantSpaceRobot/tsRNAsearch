#!/bin/bash

ncRNA_stddev=$1
tRNA_stddev=$2
gtf=$3

### Generate distribution scores
DistributionScore.R \
	$ncRNA_stddev \
	Everything_ncRNAs_cond1-vs-cond2 \
	$gtf &
DistributionScore.R \
	$tRNA_stddev \
	Everything_tRNAs_cond1-vs-cond2 \
	$gtf &
wait
### Sort the output but not the header for ncRNAs
cat Everything_ncRNAs_cond1-vs-cond2_distribution-score_all-features.tsv \
	| awk 'NR<2{print $0;next}{print $0| "sort -k15,15nr"}' \
	> Everything_ncRNAs_cond1-vs-cond2_distribution-score_all-features_sorted.tsv
### Sort the output but not the header for tRNAs
cat Everything_tRNAs_cond1-vs-cond2_distribution-score_all-features.tsv \
	| awk 'NR<2{print $0;next}{print $0| "sort -k15,15nr"}' \
	> Everything_tRNAs_cond1-vs-cond2_distribution-score_all-features_sorted.tsv
# Get features with high distribution
cat Everything_*high-distribution-score.tsv | grep -v ^feat | awk '{print $1}' \
	> High-distribution-scores_feature-names.txt
