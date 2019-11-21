# tsRNAsearch to-do list

Update after publication:
* Add batch correction: DESeq2 in-built batch correction for DESeq2 analysis. limma batch correction on RPM-normalised count data for the rest of the pipeline?
* Integrate MINTmap tRNA fragment license plates for identified tsRNAs

Update now:
* Test against Unitas, MINTmap, tDRmapper, SPORTS1.0
* Create singularity image of pipeline
* Create test to ensure Coverage-flipper script is acting correctly in human, mouse, rat
* Test setup.sh on server/virtual machine
* -R does nothing in single dataset analysis
* make a HTML report for single dataset analysis
* Flip coordinates for single dataset analysis
* In CRC_Mouse_Healthy-vs-CRC/CRC1_vs_Healthy.Results-summary HTML file, Snora68 is in results table but in none of the plot images.
* Plot Fishers test features
* Add fishers to venn diagram
* Predict which type of tsRNA each tsRNA is. e.g. 5'-tiRNA, itRNA
