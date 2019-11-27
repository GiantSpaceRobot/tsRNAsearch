# tsRNAsearch to-do list

Update after publication:
* Add batch correction: DESeq2 in-built batch correction for DESeq2 analysis. limma batch correction on RPM-normalised count data for the rest of the pipeline?
* Integrate MINTmap tRNA fragment license plates for identified tsRNAs

Update now:
* Test against Unitas, MINTmap, tDRmapper, SPORTS1.0
* Create singularity image of pipeline
* Test setup.sh on server/virtual machine
* Predict which type of tsRNA each tsRNA is. e.g. 5'-tiRNA, itRNA
* Fishers test causes problems when only 1 replicate each compared. GBM_DMSO-vs-ITE
* Remove ncRNA coverage plot from single analysis as there is text pointing to it
* Repeat analysis of CRC adenocarcinoma
