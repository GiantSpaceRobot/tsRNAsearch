#!/usr/bin/env nextflow

/*============================
Count reads using featureCount
============================*/

process FEATURE_COUNT {

    input:
    //val featureType
    file gtf
    file bamfile

    output:
    path("*.fcount"), emit: fcounts
    path("*.count"), emit: counts
   
    shell:
    '''
    featureCounts \\
        -T !task.cpus \\
        -a !gtf \\
        -o !{bamfile}.fcount \\
        !{bamfile}
    grep -v featureCounts !{bamfile}.fcount \\
        | grep -v ^Geneid \\
        | awk -v OFS="\t" "{print $1, $7}" \\
        > !{bamfile}.count
    '''
}