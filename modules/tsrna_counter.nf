#!/usr/bin/env nextflow

/*=============================
Count every tsRNA from BAM file
=============================*/

params.output_dir = 'results'
//tsRNA_out_dir = params.output_dir/"tsRNA_Data"

process TSRNA_INDIVIDUAL_COUNT {

    tag "$bamfile.simpleName"
    publishDir params.output_dir
    //publishDir tsRNA_out_dir

    input:
    file bamfile
    path(gtf)

    output:
    path("*.tsv"), emit: tsRNA_individual_counts
    //path("*.log"), emit: tsRNA_count_errors
   
    script:
    """
    samtools view -h ${bamfile} | SAM-to-tsRNAcount.py $gtf ${bamfile.simpleName}.tsv > ${bamfile.simpleName}.Errors.log
    """
}