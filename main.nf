#!/usr/bin/env nextflow
/*
 *# Author: Paul Donovan 
 *# Email: pauldonovandonegal@gmail.com
*/

nextflow.enable.dsl = 2

// Set default parameters 
params.species = 'human'
//params.skip = "no"
params.all_plots = false
//params.remove = "no"
params.layout = null
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
    nextflow run tsRNAsearch --species mouse --input_dir tsRNAsearch/ExampleData --output_dir Results

    Usage: Multi file analysis:
    nextflow run tsRNAsearch --species mouse --input_dir tsRNAsearch/ExampleData --output_dir Results

    Usage: Group comparison analysis:
    nextflow run tsRNAsearch --species mouse --input_dir tsRNAsearch/ExampleData --output_dir Results --layout tsRNAsearch/additional-files/Example_Layout.csv

    Single file analysis mandatory arguments:
    Output directory: ${params.output_dir}
    Input directory with file (/path/to/file): ${params.input_dir}

    Multi file analysis (no group comparison) mandatory arguments:
    Output directory: ${params.output_dir}
    Input directory with file (/path/to/file): ${params.input_dir}

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
version="Version:  tsRNAsearch 0.41"

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
if(!params.input_dir){
    exit 1, "Error: No input provided. Provide --input_dir to pipeline"
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
include { COUNTS_TO_COLLAPSED_COUNTS } from './modules/counts_to_collapsed_counts'
include { PREDICT_TSRNA_TYPE } from './modules/predict_tsrna_type'
include { PREDICT_TSRNA_TYPE_GROUPS } from './modules/predict_tsrna_type_groups'
include { GENERATE_RESULTS_PDF } from './modules/generate_results_pdf'
include { GENERATE_RESULTS_PDF_GROUPS } from './modules/generate_results_pdf_groups'
include { DESEQ2 } from './modules/deseq2'
include { COMBINE_DEPTHFILES } from './modules/combine_depthfiles'
include { DATA_TRANSFORMATIONS } from './modules/data_transformations'
include { DISTRIBUTION_SCORE } from './modules/distribution_score'
include { SLOPE_SCORE } from './modules/slope_score'
include { CLEAVAGE_SCORE } from './modules/cleavage_score'
include { FISHERS_METHOD } from './modules/fishers_method'
include { RESULTS_TABLE } from './modules/results_table'
include { COMBINED_SCORE } from './modules/combined_score'
include { PLOT_TSRNA_LENGTHS } from './modules/plot_tsrna_lengths'
include { PLOT_ALL_NCRNAS } from './modules/plot_all_ncrnas'
include { PLOT_TRNAS } from './modules/plot_trnas'
include { PLOT_NCRNAS } from './modules/plot_ncrnas'
include { VENN_DIAGRAM } from './modules/venn_diagram'
include { GENERATE_COUNT_DATAFRAME } from './modules/generate_count_dataframe'
include { STACKED_BARPLOTS } from './modules/stacked_barplots'
include { BARPLOTS } from './modules/barplots'
include { ORGANISE_RESULTS } from './modules/organise_results'
include { ORGANISE_RESULTS_GROUPS } from './modules/organise_results_groups'
include { PUBLISH_FILES } from './modules/publish_files'

workflow {
    main:
        // Define channels
        fastq_channel = Channel.fromPath( ["$params.input_dir/*.fastq.gz", "$params.input_dir/*.fq.gz"] ) //, "$params.input_file"] )
        PREPARE_TRNA_GTF(params.species)
        PREPARE_NCRNA_GTF(params.species)
        TRIM_READS(fastq_channel, "$params.min_read_length")
        MAKE_STAR_DB("$projectDir/DBs/${params.species}_tRNAs-and-ncRNAs-and-lookalikes.fa")  // Run process to generate DB
        STAR_ALIGN(TRIM_READS.out.trimmed_reads, MAKE_STAR_DB.out.star_index)
        BAM_COLLAPSE(STAR_ALIGN.out.bam)
        ADD_EMPTY_COUNTS(BAM_COLLAPSE.out.tRNA_almost_mapped_count, "$projectDir/additional-files/${params.species}_empty_tRNA.count")
        BAM_SPLIT(BAM_COLLAPSE.out.collapsedbam)
        FEATURE_COUNT_TRNA(BAM_SPLIT.out.bam_tRNA, PREPARE_TRNA_GTF.out.tRNA_gtf)
        FEATURE_COUNT_NCRNA(BAM_SPLIT.out.bam_ncRNA, PREPARE_NCRNA_GTF.out.ncRNA_gtf)
        SUM_COUNTS(FEATURE_COUNT_TRNA.out.counts.collect(), FEATURE_COUNT_NCRNA.out.counts.collect(), ADD_EMPTY_COUNTS.out.filled_multi_counts.collect())
        GENERATE_TRNA_DEPTH_FILES(BAM_SPLIT.out.bam_tRNA, PREPARE_TRNA_GTF.out.tRNA_gtf, PREPARE_NCRNA_GTF.out.ncRNA_gtf)
        GENERATE_MULTIMAPPER_TRNA_DEPTH_FILES(BAM_COLLAPSE.out.tRNA_almost_mapped_txt, SUM_COUNTS.out.sum_counts)
        GENERATE_NCRNA_DEPTH_FILES(BAM_SPLIT.out.bam_ncRNA, SUM_COUNTS.out.sum_counts) 
        GENERATE_DEPTHFILE_STATS(GENERATE_TRNA_DEPTH_FILES.out.depth_files.mix(GENERATE_NCRNA_DEPTH_FILES.out.depth_files))
        RAW_COUNTS_TO_PROPORTIONS(FEATURE_COUNT_TRNA.out.counts, BAM_COLLAPSE.out.tRNA_almost_mapped_count, SUM_COUNTS.out.sum_counts)
        RAW_COUNTS_TO_NORM_COUNTS(SUM_COUNTS.out.all_counts.flatten(), PREPARE_TRNA_GTF.out.tRNA_gtf, PREPARE_NCRNA_GTF.out.ncRNA_gtf)
        COUNTS_TO_COLLAPSED_COUNTS(SUM_COUNTS.out.all_counts.flatten().mix(RAW_COUNTS_TO_NORM_COUNTS.out.tpm_count))
        PREDICT_TSRNA_TYPE(GENERATE_TRNA_DEPTH_FILES.out.collect(), GENERATE_MULTIMAPPER_TRNA_DEPTH_FILES.out.collect())
        //MULTIQC(FASTQC.out.collect())

        // Plot things
        PLOT_TRNA_ALIGNMENT_LENGTH(BAM_SPLIT.out.bam_tRNA)
        PLOT_TRNA_ALL_PLOTS(GENERATE_TRNA_DEPTH_FILES.out.depth_files, PREPARE_TRNA_GTF.out.tRNA_gtf)
        GENERATE_RESULTS_PDF(PLOT_TRNA_ALIGNMENT_LENGTH.out.pdf.collect(), \
            PLOT_TRNA_ALL_PLOTS.out.pdfs.collect(), \
            SUM_COUNTS.out.sum_counts)
        if (params.all_plots){   
            PLOT_NCRNA_ALL_PLOTS(GENERATE_NCRNA_DEPTH_FILES.out.depth_files, PREPARE_NCRNA_GTF.out.ncRNA_gtf)
        }
        // END OF INDIVIDUAL FILE ANALYSIS

        // START OF GROUP COMPARISON
        if (params.layout){
            DESEQ2(COUNTS_TO_COLLAPSED_COUNTS.out.collapsed_count.collect(), "$launchDir/$params.layout", PREPARE_NCRNA_GTF.out.ncRNA_gtf)
            DATA_TRANSFORMATIONS("$launchDir/$params.layout", \
                GENERATE_TRNA_DEPTH_FILES.out.depth_files.collect(), \
                GENERATE_NCRNA_DEPTH_FILES.out.depth_files.collect(), \
                GENERATE_MULTIMAPPER_TRNA_DEPTH_FILES.out.depth_files.collect(), \
                SUM_COUNTS.out.sum_counts)
            DISTRIBUTION_SCORE(DATA_TRANSFORMATIONS.out.ncrna_stddev, DATA_TRANSFORMATIONS.out.trna_stddev, PREPARE_NCRNA_GTF.out.ncRNA_gtf)
            SLOPE_SCORE(DATA_TRANSFORMATIONS.out.depth_means, "$launchDir/$params.layout", PREPARE_NCRNA_GTF.out.ncRNA_gtf)
            CLEAVAGE_SCORE(DATA_TRANSFORMATIONS.out.depth_means, "$launchDir/$params.layout", PREPARE_NCRNA_GTF.out.ncRNA_gtf)
            FISHERS_METHOD(DATA_TRANSFORMATIONS.out.depthfiles, "$launchDir/$params.layout", PREPARE_NCRNA_GTF.out.ncRNA_gtf)
            RESULTS_TABLE(FISHERS_METHOD.out.pvalues, DISTRIBUTION_SCORE.out.combined_features, CLEAVAGE_SCORE.out.combined_features, DESEQ2.out.csvs, SLOPE_SCORE.out.combined_features, "$launchDir/$params.layout")
            COMBINED_SCORE(RESULTS_TABLE.out.tsv, "$launchDir/$params.layout")
            PLOT_TSRNA_LENGTHS(PLOT_TRNA_ALIGNMENT_LENGTH.out.txt.collect(), "$launchDir/$params.layout")
            if (params.all_plots){   
                PLOT_ALL_NCRNAS(DATA_TRANSFORMATIONS.out.ncrna_depth, "$launchDir/$params.layout", PREPARE_NCRNA_GTF.out.ncRNA_gtf)
            }
            PLOT_TRNAS(DATA_TRANSFORMATIONS.out.depthfiles, COMBINED_SCORE.out.top_features, "$launchDir/$params.layout")
            PLOT_NCRNAS(DATA_TRANSFORMATIONS.out.ncrna_depth, COMBINED_SCORE.out.top_features, "$launchDir/$params.layout", PREPARE_NCRNA_GTF.out.ncRNA_gtf)
            VENN_DIAGRAM(FISHERS_METHOD.out.top_features, DISTRIBUTION_SCORE.out.top_features, CLEAVAGE_SCORE.out.top_features, DESEQ2.out.txt, SLOPE_SCORE.out.top_features, "$launchDir/$params.layout")
            GENERATE_COUNT_DATAFRAME(COUNTS_TO_COLLAPSED_COUNTS.out.collect())
            STACKED_BARPLOTS(RAW_COUNTS_TO_PROPORTIONS.out.collect())
            BARPLOTS(GENERATE_COUNT_DATAFRAME.out, "$launchDir/$params.layout", PREPARE_NCRNA_GTF.out.ncRNA_gtf)
            PREDICT_TSRNA_TYPE_GROUPS(DATA_TRANSFORMATIONS.out.depth_means, "$launchDir/$params.layout")
            GENERATE_RESULTS_PDF_GROUPS(DESEQ2.out.pdfs, STACKED_BARPLOTS.out, BARPLOTS.out, VENN_DIAGRAM.out.pdf, COMBINED_SCORE.out.pdfs, PLOT_TRNAS.out, PLOT_NCRNAS.out, "$launchDir/$params.layout")
            ORGANISE_RESULTS_GROUPS("$launchDir/$params.output_dir", SUM_COUNTS.out.sum_counts, GENERATE_RESULTS_PDF_GROUPS.out.pdf)
            PUBLISH_FILES("$launchDir/$params.output_dir", DESEQ2.out.pdfs, STACKED_BARPLOTS.out, BARPLOTS.out, VENN_DIAGRAM.out.pdf, COMBINED_SCORE.out.pdfs, PLOT_TRNAS.out, PLOT_NCRNAS.out, "$launchDir/$params.layout", DESEQ2.out.csvs, VENN_DIAGRAM.out.tsvs, COMBINED_SCORE.out.all_txt_output)
        } else {
            ORGANISE_RESULTS("$launchDir/$params.output_dir", SUM_COUNTS.out.sum_counts, GENERATE_RESULTS_PDF.out.pdf, RAW_COUNTS_TO_PROPORTIONS.out.collect())
        }
}