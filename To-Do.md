# tsRNAsearch to-do list

Update after publication:
* Add batch correction: DESeq2 in-built batch correction for DESeq2 analysis. limma batch correction on RPM-normalised count data for the rest of the pipeline?
* Add error catching/empty file catching in python and R scripts
Update now:
* Test against Unitas, MINTmap, tDRmapper, SPORTS1.0
* Test setup.sh on server/virtual machine
* Fishers test causes problems when only 1 replicate each compared. GBM_DMSO-vs-ITE
* Correct for multiple testing for Fisher's combined p-value?  No
* Use FASTA as input?  Not right now
* Make tsRNAsearch standalone/Not require downloading stuff in the setup script
* Differentially expressed tRNAs feature plot is not ordered by top ones (or maybe it is be log2FC)
* Include txt and HTML table in output directory
* Fix output directory paths. Pipeline currently doesn't really accept absolute outdir paths
* Write tRNAs in descending score order for single file analysis
* Write pipeline stats at end of single file analysis
* Edit pipeline so it cannot remove all FASTQs on system, i.e. make rm specific
* Output simple PDF just showing QC and slope results
