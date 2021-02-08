#!/usr/bin/env nextflow
/*
 *# Author: Paul Donovan 
 *# Email: pauldonovandonegal@gmail.com
*/


nextflow.enable.dsl = 2


// Set default parameters 
params.input_file = null 
params.species = 'human'
params.skip = "no"
params.all_plots = false  // should be false
params.remove = "no"
params.layout = ""
params.input_dir = null
params.output_dir = "Results"
params.min_read_length = 16
params.version = false
params.help = false


// Print message for user
def helpMessage() {
    log.info """\
    
    ===========
    tsRNAsearch
    ===========

    Usage: Single file analysis:
    nextflow run main.nf --species mouse --input_file ExampleData/CytC_IP1.fastq.gz --output_dir Results

    Usage: Group comparison analysis:
    nextflow run main.nf --species mouse --input_dir ExampleData --output_dir Results --layout Layout.csv

    Single file analysis mandatory arguments:
    Output directory: ${params.output_dir}
    Input file: ${params.input_file}

    Group comparison analysis mandatory arguments:
    Output directory: ${params.output_dir}
    Input directory of files (/path/to/files): ${params.input_dir}
    

    Other arguments:
    Species (human/mouse/rat): ${params.species}
    Minimum read length: ${params.min_read_length}
    Print help" ${params.help}
    """
}


// Pipeline version
version="Version:  tsrna-de 0.1"


// Print message for user
def versionMessage() {
    log.info """\
    ${version}
    """
}

// Show help message and quit
if (params.help) {
    helpMessage()
    versionMessage()
    exit 0
}


// Show version and quit
if (params.version) {
    versionMessage()
    exit 0
}


// Input parameter error catching
if(!params.input_file && !params.input_dir){
    exit 1, "Error: No input provided. Provide either --input_file or --input_dir to pipeline"
}
if(params.input_dir && !params.layout){
    exit 1, "Error: No --layout file provided. See README for an example"
}
if( params.input_file && params.input_dir){
    exit 1, "Error: Conflicting inputs. Cannot supply both single FASTQ file and FASTQ input directory"
}
//if [[ "$outDir" == */ ]]; then # If outDir has a trailing slash, remove
//	outDir=$(echo "${outDir::-1}")
//fi


// Load modules (these inherit the params above if the default params are also declared in the modules)
include { PREPARE_TRNA_GTF } from './modules/prepare_trna_gtf'
include { PREPARE_NCRNA_GTF } from './modules/prepare_ncrna_gtf'
//include { FASTQC } from './modules/fastqc'
//include { MULTIQC } from './modules/multiqc'
include { MAKE_STAR_DB } from './modules/star_genome-generate'
include { TRIM_READS } from './modules/trim_galore'
include { STAR_ALIGN } from './modules/star_align'
include { BAM_COLLAPSE } from './modules/bam_collapse'
include { ADD_EMPTY_COUNTS } from './modules/add_empty_counts'
include { BAM_SPLIT } from './modules/bam_split'
include { FEATURE_COUNT_NCRNA } from './modules/count_features_ncrna'
include { FEATURE_COUNT_TRNA } from './modules/count_features_trna'
include { SUM_COUNTS } from './modules/sum_counts'
include { GENERATE_TRNA_DEPTH_FILES } from './modules/generate_trna_depth_files'
include { GENERATE_MULTIMAPPER_TRNA_DEPTH_FILES } from './modules/generate_multimapper_trna_depth_files'
include { GENERATE_NCRNA_DEPTH_FILES } from './modules/generate_ncrna_depth_files'
include { GENERATE_DEPTHFILE_STATS } from './modules/generate_depthfile_stats'
include { PLOT_TRNA_ALIGNMENT_LENGTH } from './modules/plot_trna_alignment_length'
include { PLOT_NCRNA_ALL_PLOTS } from './modules/plot_ncrna_all_plots'
include { PLOT_TRNA_ALL_PLOTS } from './modules/plot_trna_all_plots'
//include { PLOT_TRNA_SINGLE_SAMPLE_ANALYSIS } from './modules/plot_trna_single_sample_analysis'
include { PLOT_TRNA_ALL_DEPTH_PLOTS } from './modules/plot_trna_all_depth_plots'
include { RAW_COUNTS_TO_PROPORTIONS } from './modules/raw_counts_to_proportions'
include { RAW_COUNTS_TO_NORM_COUNTS } from './modules/raw_counts_to_norm_counts'
include { RAW_COUNTS_TO_COLLAPSED_COUNTS } from './modules/raw_counts_to_collapsed_counts'
include { PREDICT_TSRNA_TYPE } from './modules/predict_tsrna_type'
include { ORGANISE_RESULTS } from './modules/organise_results'
include { GENERATE_RESULTS_PDF } from './modules/generate_results_pdf'
include { DESEQ2 } from './modules/deseq2'

//include { TSRNA_INDIVIDUAL_COUNT } from './modules/tsrna_counter'
//include { TSRNA_DESEQ } from './modules/tsrna_deseq2'


workflow {
    main:
        // Define channels
        //ncRNA_gtf = Channel.fromPath("$projectDir/DBs/${params.species}_ncRNAs_relative_cdhit.gtf")
        fastq_channel = Channel.fromPath( ["$params.input_dir/*.fastq.gz", "$params.input_dir/*.fq.gz"] ) //, "$params.input_file"] )
        //fastq_channel.view()
        //exit 1
        PREPARE_TRNA_GTF(params.species)
        PREPARE_NCRNA_GTF(params.species)
        //PREPARE_TRNA_GTF.out.tRNA_gtf.view()
        //FASTQC(fastq_channel)
        //MULTIQC(FASTQC.out.collect())
        TRIM_READS(fastq_channel, "$params.min_read_length")
        //TRIM_READS.out.trimmed_reads.view()
        MAKE_STAR_DB("$projectDir/DBs/${params.species}_tRNAs-and-ncRNAs_relative_cdhit.fa")  // Run process to generate DB
        STAR_ALIGN(TRIM_READS.out.trimmed_reads, MAKE_STAR_DB.out.star_index)
        //PREPARE_NCRNA_GTF.out.ncRNA_gtf.view()
        BAM_COLLAPSE(STAR_ALIGN.out.bam)
        ADD_EMPTY_COUNTS(BAM_COLLAPSE.out.tRNA_almost_mapped_count, "$projectDir/additional-files/${params.species}_empty_tRNA.count")
        BAM_SPLIT(BAM_COLLAPSE.out.collapsedbam)
        FEATURE_COUNT_TRNA(BAM_SPLIT.out.bam_tRNA, PREPARE_TRNA_GTF.out.tRNA_gtf)
        FEATURE_COUNT_NCRNA(BAM_SPLIT.out.bam_ncRNA, PREPARE_NCRNA_GTF.out.ncRNA_gtf)
        SUM_COUNTS(FEATURE_COUNT_TRNA.out.counts.collect(), FEATURE_COUNT_NCRNA.out.counts.collect(), ADD_EMPTY_COUNTS.out.filled_multi_counts.collect())
        GENERATE_TRNA_DEPTH_FILES(BAM_SPLIT.out.bam_tRNA, SUM_COUNTS.out.sum_counts)
        GENERATE_MULTIMAPPER_TRNA_DEPTH_FILES(BAM_COLLAPSE.out.tRNA_almost_mapped_txt, SUM_COUNTS.out.sum_counts)
        GENERATE_NCRNA_DEPTH_FILES(BAM_SPLIT.out.bam_ncRNA, SUM_COUNTS.out.sum_counts) 
        GENERATE_DEPTHFILE_STATS(GENERATE_TRNA_DEPTH_FILES.out.depth_files.mix(GENERATE_NCRNA_DEPTH_FILES.out.depth_files))
        RAW_COUNTS_TO_PROPORTIONS(FEATURE_COUNT_TRNA.out.counts, BAM_COLLAPSE.out.tRNA_almost_mapped_count, SUM_COUNTS.out.sum_counts)
        RAW_COUNTS_TO_NORM_COUNTS(SUM_COUNTS.out.all_counts.flatten(), SUM_COUNTS.out.sum_counts)
        RAW_COUNTS_TO_COLLAPSED_COUNTS(SUM_COUNTS.out.all_counts.flatten().mix(RAW_COUNTS_TO_NORM_COUNTS.out.rpm_count))
        PREDICT_TSRNA_TYPE(GENERATE_TRNA_DEPTH_FILES.out.collect(), GENERATE_MULTIMAPPER_TRNA_DEPTH_FILES.out.collect())

        // Plot things
        PLOT_TRNA_ALIGNMENT_LENGTH(BAM_SPLIT.out.bam_tRNA)
        PLOT_TRNA_ALL_PLOTS(GENERATE_TRNA_DEPTH_FILES.out.depth_files, PREPARE_TRNA_GTF.out.tRNA_gtf)
        if (params.all_plots){   
            //PLOT_NCRNA_ALL_PLOTS(GENERATE_NCRNA_DEPTH_FILES.out.depth_files, PREPARE_NCRNA_GTF.out.ncRNA_gtf)
        }
        GENERATE_RESULTS_PDF(PLOT_TRNA_ALIGNMENT_LENGTH.out.pdf.collect(), PLOT_TRNA_ALL_PLOTS.out.pdfs.collect(), SUM_COUNTS.out.sum_counts)
        // END OF INDIVIDUAL FILE ANALYSIS

        // START OF GROUP COMPARISON
        DESEQ2(RAW_COUNTS_TO_COLLAPSED_COUNTS.out.collapsed_count.collect(), "$launchDir/$params.layout", PREPARE_NCRNA_GTF.out.ncRNA_gtf)
        // Organise results directory
        // This is not waiting until all processes are finished. I need to explicitly give it output of processes
        //ORGANISE_RESULTS("$launchDir/$params.output_dir", SUM_COUNTS.out.sum_counts, PLOT_TRNA_ALL_PLOTS.out.pdfs)

        //TSRNA_INDIVIDUAL_COUNT(SAM_SPLIT_AND_SAM2BAM.out.bam_tRNA, PREPARE_TRNA_GTF.out.tRNA_gtf)
        //PREPARE_TRNA_GTF.out.tRNA_gtf.view()
        //TSRNA_DESEQ("$launchDir/${params.layout}", TSRNA_INDIVIDUAL_COUNT.out.tsRNA_individual_counts.collect())



}

