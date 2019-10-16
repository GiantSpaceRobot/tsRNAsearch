# tsRNAsearch to-do list

Update after publication:
* Add batch correction: DESeq2 in-built batch correction for DESeq2 analysis. limma batch correction on RPM-normalised count data for the rest of the pipeline?
* Should I normalise my Fisher method p-values by the length of each feature? Or the number of nucleotides that had more than 0 coverage?

Update now:
* Integrate MINTmap tRNA fragment license plates for identified tsRNAs
* Use GRCh38 ncRNA sequences
* Place the pipeline in a docker container?
* Create one .txt or .csv with all summary results.
	e.g. ValCAC	DESeq2 Log2FC	DESeq2 padj	Distribution score	Cleavage score	Overall fisher p-val
* Remove genomes from setup
* Change Log2FC and padj cut-offs for DESeq2? If so, also edit Volcano plot cut-offs
