#!/bin/bash

# Setup for tiRNA pipeline

python 2.7.15rc1
	numpy
	sys
Rscript 3.4.4
	ggplot2
	gplots

samtools 1.7 # Requires htslib 1.7-2
trim_galore 0.4.4_dev
fastqc v0.11.5
hisat2 2.1.0
featureCounts 1.6.3




