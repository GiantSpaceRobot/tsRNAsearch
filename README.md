# tsRNAsearch

A Nextflow DSL2 pipeline for the identification, quantification and analysis of ncRNAs (especially tRNA fragments) in small/miRNA-seq datasets

Many ncRNA identification pipelines are based on differential gene expression. This is a useful method but can fall short if portions of the differentially expressed ncRNAs are degraded/not present, thus reducing the differential expression identification power. tsRNAsearch addresses this issue by applying a standard differential expression analysis and four additional methods to identify fragmented ncRNAs, especially tRNA fragments.

## Quickstart
### Installation (via new conda environment)
Run the code shown in the window below. This will create a conda environment and install the required tools inside the environment. 

NOTE: please use conda 4.6 or greater

```
git clone https://github.com/GiantSpaceRobot/tsRNAsearch.git
conda env create -f tsRNAsearch/environment.yml
conda activate tsrnasearch_env
Rscript tsRNAsearch/bin/InstallPackages.R
```

### Running tsRNAsearch
We have supplied data to ensure that the pipeline is functioning correctly.
To run the following, make sure you are outside the tsRNAsearch directory.

Usage: Single file analysis:
```
nextflow run tsRNAsearch --species mouse --input_dir tsRNAsearch/ExampleData --output_dir Results
```

Usage: Multi file analysis:
```
nextflow run tsRNAsearch --species mouse --input_dir tsRNAsearch/ExampleData --output_dir Results
```
On finishing either the single or multi-file analysis, the pipeline will populate the Results directory with results files corresponding to each sample.

Usage: Group comparison analysis:
```
nextflow run tsRNAsearch --species mouse --input_dir tsRNAsearch/ExampleData --output_dir Results --layout tsRNAsearch/additional-files/Example_Layout.csv
```
On finishing the group comparison analysis, the pipeline will populate the Results directory with results files corresponding to each sample, a PDF summarising the comparison, and a directory containing comparison files (e.g. DESeq2 results). [Example of pipeline results](https://giantspacerobot.github.io/tsRNAsearch_ExampleOutput/)

## More Information
#### ncRNA identification methods:
* Slope algorithm
* Fisher's method for combining p-values
* Distribution algorithm
* Cleavage algorithm
* DESeq2

### tsRNAsearch Steps and Parameters 
#### Steps:
* Execute standard small RNA-seq steps on each dataset
  1. trim\_galore (cutadapt and fastqc) 
    * Adapter removal, read trimming, and quality check
  2. STAR
    * Read alignment to a non-coding RNA database
  3. FeatureCounts
    * Count raw reads
  4. Data processing
* Data transformations/reformatting
* Run SlopeScore.R
* Run Pvalue_generator.R
* Run DistributionScore.R
* Run CleavageScore.R
* Run DESeq2
* Generate reports (PDFs, CSVs, and text reports)
#### Parameters:
* -h *Print the usage and options information*
* --species *Analyse data against '__human__', '__mouse__', or '__rat__'? {default: __human__}*
* --input_dir *Directory containing the files for analysis (file formats: FASTQ or gzipped FASTQ). Directory should have no other contents.*
* --output_dir *Output directory for the results files*
* --layout *CSV file containing file names and file groups (see example in additional-files/)*
* --min_read_length *Minimum read length (default: __16__ bp)*
* --all_plots *Plot all features? __true__/__false__ {default: __false__}*
* --help *Print help and exit*
* --version *Print version and exit*

### Methods to Address Difficulty in Correctly Identifying tRNA reads
There are 625 tRNA genes in the human genome. Each tRNA species has multiple genes coding for it, many of which are identical or vastly similar (e.g. Arginine). Therefore, correctly matching a read with its tRNA gene of origin is difficult. Two methods have been implemented to address this:
* CD-HIT was used to extract representative sequences from all ncRNA sequences. The ncRNA STAR database was built using these representative sequences.
* The STAR output (SAM file) is processed in such a way that all reads mapping to multiple database sequences from the same ncRNA species (e.g. read maps to multiple Proline-CCC genes) are collapsed into a single SAM entry. This single SAM entry is thereby considered a single-mapping read and is not removed by downstream processing steps.

### Adding additional species
You can add any species you want with relatively little effort. Click [here](https://github.com/GiantSpaceRobot/tsRNAsearch_add-new-species) for more information.

### More info
You can test if the R packages are installed correctly in the conda environment using:
```
# You must be in the conda env to run this
Rscript tsRNAsearch/bin/TestPackages.R
```
This will produce errors if all packages are not found. You may need to manually install the problem packages if so.

## Contributors
* Paul Donovan, PhD
* Natalie McHale

## License
This project is licensed under the MIT License.

