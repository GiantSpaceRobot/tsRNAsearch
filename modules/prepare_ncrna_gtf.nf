#!/usr/bin/env nextflow

/*====================
Prepare ncRNA GTF file
Place it in a channel
====================*/

process PREPARE_NCRNA_GTF {

    //tag 

    input:
    val species

    output:
    path("*ncRNA*.gtf"), emit: ncRNA_gtf
   
    script:
    """
    cp "$projectDir/DBs/${species}_ncRNAs_relative_cdhit.gtf" .
    """
}