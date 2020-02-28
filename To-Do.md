# tsRNAsearch to-do list

Update after publication:
* Add batch correction: DESeq2 in-built batch correction for DESeq2 analysis. limma batch correction on RPM-normalised count data for the rest of the pipeline?

Update now:
* Test against Unitas, MINTmap, tDRmapper, SPORTS1.0
* Test setup.sh on server/virtual machine
* Fishers test causes problems when only 1 replicate each compared. GBM_DMSO-vs-ITE
* Add error catching/empty file catching in python and R scripts
* Correct for multiple testing for Fisher's combined p-value?
* Make plotting scripts faster?
* Use FASTA as input?
* Replace full path to Plots with relative path in HTML, allowing HTML to move between computers as long as Plots dir is there 
* change min read length to 16
* Fix HTMLs
* Sort gene lists prior to plotting in Bedgraph_plotter_DEGs.R script so that plots are sorted by rank, not alphabetically
