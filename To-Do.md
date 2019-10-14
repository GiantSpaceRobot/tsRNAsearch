# tsRNAsearch to-do list

Update after publication:
* Add batch correction: DESeq2 in-built batch correction for DESeq2 analysis. limma batch correction on RPM-normalised count data for the rest of the pipeline?
* Should I normalise my Fisher method p-values by the length of each feature? Or the number of nucleotides that had more than 0 coverage?

Update now:
* Integrate MINTmap tRNA fragment license plates for identified tsRNAs
* Use GRCh38 ncRNA sequences
* Place the pipeline in a docker container?
* Count total number of tsRNA/miRNA reads as a proportion of total reads in dataset. Make a pie chart (with relative size? i.e. pie chart should be bigger if more reads mapped in total in this condition) or bar plot with relative amounts. Compare disease vs control
* Distribution and cleavage algorithms need to pick up transcription in one condition but not the other as only p-value method picks these up
* Create volcano plot with DESeq2 script
* Create one .txt or .csv with all summary results.
	e.g. ValCAC	DESeq2 Log2FC	DESeq2 padj	Distribution score	Cleavage score	Overall fisher p-val
* Remove genomes from setup
