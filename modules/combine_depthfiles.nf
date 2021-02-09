#!/usr/bin/env nextflow

/*=================
Combine depthfiless
=================*/


process COMBINE_DEPTHFILES {

    //tag "$trna_depthfile.simpleName"

    input:
    path trna_depthfile
    path ncrna_depthfile
    path multi_trna_depthfile
    path mapped_read_counts

    output:
    path("*tsv"), emit: tsvs

    script:
    """
    cat ${mapped_read_counts} | grep -v ^Sample | while read line
    do
        my_sample=\$(echo \$line | awk -F '_trimmed' '{print \$1}')
        echo Combining files for \$my_sample
    done
    touch my.tsv
    """
}