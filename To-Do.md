# tsRNAsearch to-do list

* Add batch correction: DESeq2 in-built batch correction for DESeq2 analysis. limma batch correction on RPM-normalised count data for the rest of the pipeline?
* Generate p-value for features in cleavage and distribution algorithms. Check paper draft for p-value idea. Wilcoxon signed-rank test?
* Integrate MINTmap tRNA fragment license plates for identified tsRNAs
* Use GRCh38 ncRNA sequences
* Place the pipeline in a docker container?
* Make Distribution algorithm slightly more stringent
* Count total number of tsRNA/miRNA reads as a proportion of total reads in dataset. Make a pie chart (with relative size? i.e. pie chart should be bigger if more reads mapped in total in this condition) or bar plot with relative amounts. Compare disease vs control
* Set fastp max threads to 16
* Distribution and cleavage algorithms need to pick up transcription in one condition but not the other as only p-value method picks these up
* Should I normalise my Fisher method p-values by the length of each feature? Or the number of nucleotides that had more than 0 coverage?
* Is % area difference alone a good score?
	Probably not because it weights tiny transcription and huge transcription levels the same. Hence the distribution score.
* P-value calculation: I need to implement something to deal with single replicate comparisons.  
* Get mean rpm for each condition for each feature, plot condition 1 vs condition 2
	This is really just a DE experiment. Same can be achieved with volcano plot
