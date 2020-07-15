# tsRNAsearch to-do list

Update after publication:
* Add batch correction: DESeq2 in-built batch correction for DESeq2 analysis. limma batch correction on RPM-normalised count data for the rest of the pipeline?
* Add error catching/empty file catching in python and R scripts
* Use FASTA as input?
Update now:
* Fix output directory paths. Pipeline currently doesn't really accept absolute outdir paths
* Edit pipeline so it cannot remove all FASTQs on system, i.e. make rm specific
* Shorten file names, they are too long. This causes an error when HTMLs render long PDF paths
* Add multimapper counts to main count file so DESeq2 will use it. I'll need to add all tRNA species, even if they have 0 reads mapping as the number of features and order are vitally important. 
* Make script to analyse tsRNA-aligned.sam and output percentage of reads that are 5'-tiRNA, 5'-tRF, etc. Generate standard deviation for this and compare both conditions using t-test. Show this in boxplot 
* Insert t-test table into HTML and PDF (if possible)
* Bug in Barplots.R: Wrong label randomly being assigned to individual barplots PDF
* Remove markdown Rscript from this dir
* Only top 4 tsRNAs shown in combined score FeaturePlot
* Get rid of long and short PDF. Just the long one is fine.
* Summary report table in HTML doesn't show slope
* Remove the 'More Results' section of HTML
* Combined score not showing in HTML
* Correct Fisher's combined method for multiple testing?
* Should SlopeScore.R have penalty function? It currently does not
