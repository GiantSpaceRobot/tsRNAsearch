# tsRNAsearch to-do list

Update after publication:
* Add batch correction: DESeq2 in-built batch correction for DESeq2 analysis. limma batch correction on RPM-normalised count data for the rest of the pipeline?
* Add error catching/empty file catching in python and R scripts
* Use FASTA as input?
Update now:
* Test setup.sh on server/virtual machine
* Make tsRNAsearch standalone/Not require downloading stuff in the setup script
* Fix output directory paths. Pipeline currently doesn't really accept absolute outdir paths
* Edit pipeline so it cannot remove all FASTQs on system, i.e. make rm specific
* Shorten file names, they are too long
* Add multimapper counts to main count file so DESeq2 will use it. I'll need to add all tRNA species, even if they have 0 reads mapping as the number of features and order are vitally important. 
* Make script to analyse tsRNA-aligned.sam and output percentage of reads that are 5'-tiRNA, 5'-tRF, etc. Generate standard deviation for this and compare both conditions using t-test. Show this in boxplot 
* Insert t-test table into HTML and PDF (if possible)
