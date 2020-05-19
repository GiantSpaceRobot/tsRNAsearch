# tsRNAsearch to-do list

Update after publication:
* Add batch correction: DESeq2 in-built batch correction for DESeq2 analysis. limma batch correction on RPM-normalised count data for the rest of the pipeline?
* Add error catching/empty file catching in python and R scripts
* Use FASTA as input?
Update now:
* Test against Unitas, MINTmap, tDRmapper, SPORTS1.0
* Test setup.sh on server/virtual machine
* Fishers test causes problems when only 1 replicate each compared. GBM_DMSO-vs-ITE
* Make tsRNAsearch standalone/Not require downloading stuff in the setup script
* Differentially expressed tRNAs feature plot is not ordered by top ones (or maybe it is be log2FC)
* Include txt and HTML table in output directory
* Fix output directory paths. Pipeline currently doesn't really accept absolute outdir paths
* Write tRNAs in descending score order for single file analysis
* Write pipeline stats at end of single file analysis
* Edit pipeline so it cannot remove all FASTQs on system, i.e. make rm specific
* Shorten file names, they are too long
* Make an excel file with two columns, one for each condition, with RPM count for every feature.
* Add barplot breakdown of ncRNAs for each sample
* Add header to final FCount RPM count file 
* Add multimapper counts to main count file so DESeq2 will use it. I'll need to add all tRNA species, even if they have 0 reads mapping as the number of features and order are vitally important. 
* Fancy HTML generation is broken 
