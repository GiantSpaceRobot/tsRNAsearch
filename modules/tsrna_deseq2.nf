#!/usr/bin/env nextflow

/*========================
Carry out DESeq2 on tsRNAs
========================*/

params.output_dir = 'results'
//tsRNA_out_dir = params.output_dir/"tsRNA_Data"

process TSRNA_DESEQ {
    conda '/home/paul/Documents/Applications/Miniconda2/miniconda2/envs/tsrnasearch_env'
    //tag "$tsvfiles"
    publishDir params.output_dir
    //publishDir tsRNA_out_dir

    input:
    path layout
    file tsvfiles

    output:
    path("*.csv"), emit: tsRNA_de_results
    path("*.pdf"), emit: tsRNA_pdfs
   
    script:
    """
    conda env list
    tsRNA-Conditional-Differences.R ${layout} ${tsvfiles} > DESeq2-Errors.log
    """
}