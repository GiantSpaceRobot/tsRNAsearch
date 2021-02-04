#!/usr/bin/env nextflow

/*===================
Prepare GTF file
Place it in a channel
===================*/

process PREPARE_TRNA_GTF {

    //tag 

    input:
    val species

    output:
    path("*tRNA*.gtf"), emit: tRNA_gtf
   
    script:
    """
    cp "$projectDir/DBs/${species}_tRNAs_relative_cdhit.gtf" .
    """
}