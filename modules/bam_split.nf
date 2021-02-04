#!/usr/bin/env nextflow

/*============================
Split BAM file by feature type
============================*/

process BAM_SPLIT {

    input:
    path bamfile

    output:
    path("*ncRNAs.bam"), emit: bam_ncRNA
    path("*tRNAs.bam"), emit: bam_tRNA
    //path("*.bam"), emit: bams

    script:
    """
    BAM-split.sh ${bamfile} > BAM_split.log
    """
}