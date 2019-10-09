# tsRNAsearch to-do list

* Add batch correction: DESeq2 in-built batch correction for DESeq2 analysis. limma batch correction on RPM-normalised count data for the rest of the pipeline?
* Generate p-value for features in cleavage and distribution algorithms. Check paper draft for p-value idea. Wilcoxon signed-rank test?
* Integrate MINTmap tRNA fragment license plates for identified tsRNAs
* Use GRCh38 ncRNA sequences
* Place the pipeline in a docker container?
* Make Distribution algorithm slightly more stringent
* Make cleavage score more stringent (it always identifies more features than the other two methods)
* Count total number of tsRNA/miRNA reads as a proportion of total reads in dataset. Make a pie chart. Compare disease vs control
