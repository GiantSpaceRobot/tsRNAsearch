# tsRNAsearch

A pipeline for the identification, quantification and analysis of ncRNAs (especially tRNA fragments) in small/miRNA-seq datasets

Many ncRNA identification pipelines are based on differential gene expression. This is a useful method but can fall short if portions of the differentially expressed ncRNAs are degraded/not present, thus reducing the differential expression identification power. tsRNAsearch addresses this issue by applying a standard differential expression analysis and two additional methods to identify fragmented ncRNAs, especially tRNA fragments.

## Quickstart
There are two options for installation:
* Via conda environment
* Via global system install (requires sudo/root access)

### Installation (via new conda environment)
Run the code shown in the window below. This will create a conda environment and install the required tools inside the environment. 

NOTE: please use conda 4.6 or greater

```
git clone https://github.com/GiantSpaceRobot/tsRNAsearch.git
cd tsRNAsearch
bash setup_conda.sh
conda activate tsrnasearch_env
#Please indicate the species you wish to analyse using the -s parameter.
#Options are human/mouse/rat/all.
bash bin/setup_conda_2.sh -s all
```

### Installation (via root/sudo)
Alternatively, run the setup.sh script to install required software systemwide:

```
git clone https://github.com/GiantSpaceRobot/tsRNAsearch.git
cd tsRNAsearch
#Please indicate the species you wish to analyse using the -s parameter.
#Options are human/mouse/rat/all.
sudo bash setup_global.sh -s all
```

This script will install the following pieces of software at root level:

* fastqc
* cutadapt
* python/pip (and python modules)
* R/Rscript (and R libraries)
* samtools

### Running tsRNAsearch

tsRNAsearch can now be executed from the tsRNAsearch directory.
To run it from anywhere in your system, please run the following to create absolute paths in tsRNAsearch and add the program to your path:

```
#Create absolute path names in tsRNAsearch
bash bin/SetPath.sh

#Add ~/bin to your PATH to allow your system to find tsRNAsearch
#To do so, add the following line at the end of your .bashrc
export PATH="$PATH:~/bin"

#Copy or link tsRNAsearch to ~/bin
mkdir -p ~/bin # If ~/bin does not exist, create it
cd ~/bin
ln -s /path/to/tsRNAsearch/tsRNAsearch ./ # Create the symbolic link
```

### Analysing a single RNA-seq dataset
We have supplied data to test that the pipeline is functioning correctly:

```
bash tsRNAsearch -s mouse -f ExampleData/CytC_IP1.fastq.gz -o CytC_Results > Run-report.log
```

On finishing the run, the pipeline will produce a run report as direct output (Run-report.log). Results files will appear in *CytC\_Results/Data\_and\_Plots/*.

### Comparing two conditions (e.g. control vs treatment)
```
bash tsRNAsearch -s mouse -d ExampleData/ -e additional-files/Example_ExperimentLayout.csv -o MyResults > Run-report.log 
```

On finishing the run, PDF and HTML reports will appear in the *MyResults* directory. In addition, *MyResults* will contain directories named *Data* and *Plots* containing data and plots. [Example of pipeline results](https://giantspacerobot.github.io/tsRNAsearch_ExampleOutput/)

## More Information
#### ncRNA identification methods:
* Slope algorithm
* Fisher's combined p-value method
* Distribution algorithm
* Cleavage algotithm
* DESeq2

### tsRNAsearch Steps and Parameters 
#### Steps:
* Execute tsRNAsearch.sh on each dataset
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
* Generate reports (HTML, PDFs, CSVs, and text reports)
#### Parameters:
* -h *Print the usage and options information*
* -s *Analyse data against '__human__', '__mouse__', or '__rat__'? {default: __human__}*
* -d *Directory containing the files for analysis (file formats: FASTQ or gzipped FASTQ). Directory should have no other contents.*
* -f *Single-end file for analysis*
* -o *Output directory for the results files*
* -e *CSV file containing file names and file groups (see example in additional-files/)*
* -t *Number of threads to use {default is to calculate the number of processors and use 75%}*
* -l *Minimum read length (default: __16__ bp)*
* -A *Plot all features? __yes__/__no__ {default: __no__}*
* -S *Skip pre-processing of data (i.e. skip trim_galore) {default: __no__}*
* -R *Remove all unnecessary/intermediate files {default: __no__}*

### Methods to Address Difficulty in Correctly Identifying tRNA reads
There are 625 tRNA genes in the human genome. Each tRNA species has multiple genes coding for it, many of which are identical or vastly similar (e.g. Arginine). Therefore, correctly matching a read with its tRNA gene of origin is difficult. Two methods have been implemented to address this:
* CD-HIT was used to extract representative sequences from all ncRNA sequences. The ncRNA STAR database was built using these representative sequences.
* The STAR output (SAM file) is processed in such a way that all reads mapping to multiple database sequences from the same ncRNA species (e.g. read maps to multiple Proline-CCC genes) are collapsed into a single SAM entry. This single SAM entry is thereby considered a single-mapping read and is not removed by downstream processing steps.

### Running tsRNAsearch using FASTA as input instead of FASTQ
To use the pipeline with FASTA data, you need to skip preprocessing by using the '-S yes' parameter. 

### Adding additional species
You can add any species you want with relatively little effort. Click [here](https://github.com/GiantSpaceRobot/tsRNAsearch_add-new-species) for more information.

## Contributors
* Paul Donovan, PhD
* Natalie McHale

## License
This project is licensed under the MIT License.

