# tsRNAsearch to-do list

* Add batch correction: DESeq2 in-built batch correction for DESeq2 analysis. limma batch correction on RPM-normalised count data for the rest of the pipeline?
* Add multimapper counts to main count file so DESeq2 will use it. I'll need to add all tRNA species, even if they have 0 reads mapping as the number of features and order are vitally important. 
* Insert t-test table into HTML and PDF (if possible)
* Remove markdown Rscript from this dir
* Only top 4 tsRNAs shown in combined score FeaturePlot
* Summary report table in HTML doesn't show slope
* Combined score not showing in HTML: update ResultsTable.R to include Combined scores
* Correct Fisher's combined method for multiple testing?
* DESeq2 tRNA plots not sorted correctly (mouse mito tRNA pipeline version) (quick)
* Add mito tRNAs to human and rat DBs
* Combined score plots fail to generate when n=1
