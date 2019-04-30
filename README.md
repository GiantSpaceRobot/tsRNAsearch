# tsRNAsearch

A pipeline for the identification, quantification and analysis of ncRNAs (especially tRNA fragments) in small/miRNA-seq datasets

Many ncRNA identification pipelines are based on differential gene expression. This is a useful method but can fall short if portions of the differentially expressed ncRNAs are degraded, thus reducing the DE identification power. tsRNAsearch addresses this issue by applying a standard DE analysis and two separate (but related) methods to identify fragmented ncRNAs, especially tsRNAs.

## INSTALL
```
#chmod 755 setup.sh
sudo ./setup.sh -g human # (human/mouse/both)
```
Add tsRNAsearch.sh and tsRNAsearch\_DE.sh to your path
## Quickstart
### Analysing a single RNA-seq dataset
```
tsRNAsearch.sh -g human -s ExampleData/CytC_IP1.fastq.gz -o CytC_Results -t 1
```
### Comparing two conditions (e.g. control vs treatment)
```
tsRNAsearch_DE.sh -g human -d ExampleData/ -e additional-files/Example_ExperimentLayout.csv -o MyResults -t 1 
```
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
* HISAT2 (Reads mapped to ncRNA database. Unmapped reads mapped to whole genome)
* FeatureCounts
* Data processing
#### Parameters:
* -h *Print the usage and options information*
* -g *Analyse data against 'human' or 'mouse'? {default: human}*
* -s *Single-end file for analysis*
* -o *Output directory for the results and log files*
* -A *Plot all features? yes/no {default: yes}*
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
* -g *Analyse datasets against 'human' or 'mouse'? {default: human}*
* -d *Directory containing the files for analysis (file formats: FASTQ or gzipped FASTQ). Directory should have no other contents.*
* -o *Output directory for the results and log files*
* -e *Optional (but recommended) CSV file containing file names and file groups (see examples in additional-files/)*
* -t *Number of threads to use {default is to calculate the number of processors and use 75%}*
* -A *Plot all features? yes/no {default: no (only plot differentially expressed features)}*

## Contributors
* Paul Donovan, PhD

## License
This project is licensed under the MIT License.

