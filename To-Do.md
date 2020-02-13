# tsRNAsearch to-do list

Update after publication:
* Add batch correction: DESeq2 in-built batch correction for DESeq2 analysis. limma batch correction on RPM-normalised count data for the rest of the pipeline?
* Integrate MINTmap tRNA fragment license plates for identified tsRNAs

Update now:
* Test against Unitas, MINTmap, tDRmapper, SPORTS1.0
* Test setup.sh on server/virtual machine
* Fishers test causes problems when only 1 replicate each compared. GBM_DMSO-vs-ITE
* Add error catching/empty file catching in python and R scripts
* Correct for multiple testing for Fisher's combined p-value?
* Bug in tRNA_Alignment_Length_multi-replicates.R: Can't bind because some arguments have the same name. e.g. Batch3_PD-Dementia_vs_PD-None
