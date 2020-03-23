# tsRNAsearch to-do list

Update after publication:
* Add batch correction: DESeq2 in-built batch correction for DESeq2 analysis. limma batch correction on RPM-normalised count data for the rest of the pipeline?

Update now:
* Test against Unitas, MINTmap, tDRmapper, SPORTS1.0
* Test setup.sh on server/virtual machine
* Fishers test causes problems when only 1 replicate each compared. GBM_DMSO-vs-ITE
* Add error catching/empty file catching in python and R scripts
* Correct for multiple testing for Fisher's combined p-value?
* Use FASTA as input?
* Output or convert venn diagram to PDF and add to summary PDF
* DE tRNAs are in slightly wrong order in Test 4
* Integrate slope method into single sample analysis
* Make tsRNAsearch standalone/Not require downloading stuff in the setup script
* Run Cleavage, Distribution, Slope methods on multi-mappers too
