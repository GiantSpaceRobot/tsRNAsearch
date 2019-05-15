# tsRNAsearch

A pipeline for the identification, quantification and analysis of ncRNAs (especially tRNA fragments) in small/miRNA-seq datasets

Many ncRNA identification pipelines are based on differential gene expression. This is a useful method but can fall short if portions of the differentially expressed ncRNAs are degraded, thus reducing the DE identification power. tsRNAsearch addresses this issue by applying a standard DE analysis and two separate (but related) methods to identify fragmented ncRNAs, especially tsRNAs.

## Quickstart
### Installation
Run the setup.sh script to install required software:

```
#Please indicate the genome you wish to analyse using the -g parameter. Options are human/mouse/both.
sudo ./setup.sh -g both
```

This script will install the following pieces of software at root level:

* FastQC
* Cutadapt
* STAR
* python/pip (and python modules)
* R/Rscript (and R libraries)
* samtools

Add tsRNAsearch.sh and tsRNAsearch\_DE.sh to your path

We have supplied data to test that the pipeline is functioning correctly:

### Analysing a single RNA-seq dataset
```
tsRNAsearch.sh -g mouse -s ExampleData/CytC_IP1.fastq.gz -o CytC_Results -t 1
```

On finishing the run, search for files in *CytC\_Results/Data\_and\_Plots/* to confirm completion of tsRNAsearch run. 

### Comparing two conditions (e.g. control vs treatment)
```
tsRNAsearch_DE.sh -g mouse -d ExampleData/ -e additional-files/Example_ExperimentLayout.csv -o MyResults -t 1 
```

On finishing the run, search for files in *MyResults/Plots/* to confirm completion of tsRNAsearch run.

## More Information
#### ncRNA identification methods:
* DESeq2
* Distribution score
* Cleavage score

tsRNAsearch\_DE.sh executes tsRNAsearch.sh on each RNA-seq file, concatenates the results and uses the three aformentioned methods to identify ncRNAs (both fragmented and full), and genes.

### tsRNAsearch.sh 
#### Steps:
* Trim\_galore 
* FastQC
* STAR (Reads mapped to ncRNA database. Unmapped reads mapped to whole genome)
* FeatureCounts
* Data processing
#### Parameters:
* -h *Print the usage and options information*
* -g *Analyse data against '__human__' or '__mouse__'? {default: __human__}*
* -s *Single-end file for analysis*
* -o *Output directory for the results and log files*
* -A *Plot all features? __yes__/__no__ {default: __yes__}*
* -t *Number of threads to use {default is to calculate the number of processors and use 75%}*

### tsRNAsearch\_DE.sh 
#### Steps:
* Execute tsRNAsearch.sh on each dataset
* Data transformations/reformatting
* Run DESeq2
* Run DistributionScore.R
* Run CleavageScore.R
* Generate PDFs, CSVs and text reports
#### Parameters
* -h *Print the usage and options information*
* -g *Analyse datasets against '__human__' or '__mouse__'? {default: __human__}*
* -d *Directory containing the files for analysis (file formats: FASTQ or gzipped FASTQ). Directory should have no other contents.*
* -o *Output directory for the results and log files*
* -e *Optional (but recommended) CSV file containing file names and file groups (see examples in additional-files/)*
* -t *Number of threads to use {default is to calculate the number of processors and use 75%}*
* -A *Plot all features? __yes__/__no__ {default: __no__ (only plot differentially expressed features)}*

### Methods to Address Difficulty in Correctly Identifying tRNA reads
There are 625 tRNA genes in the human genome. Each tRNA species (e.g. Arginine) has multiple genes coding for it, many of which are identical or vastly similar. Therefore, correctly matching a read with its tRNA gene of origin is difficult. Two methods have been implemented to address this:
* CD-HIT was used to extract representative sequences from all ncRNA sequences. The ncRNA STAR database was built using these representative sequences.
* The STAR output (SAM file) is processed in such a way that all reads mapping to multiple database sequences from the same ncRNA species (e.g. read maps to multiple Proline-CCC genes) are collapsed into a single SAM entry. This single SAM entry is thereby considered a single-mapping read and is not removed by downstream processing steps.

## Contributors
* Paul Donovan, PhD

## License
This project is licensed under the MIT License.

