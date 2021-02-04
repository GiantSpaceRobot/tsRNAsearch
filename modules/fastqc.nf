#!/usr/bin/env nextflow

/*========
Run FastQC
========*/

process FASTQC {

    input:
    file reads

    output:
    file "*_fastqc.{zip,html}"

    script:
    """
    fastqc $reads
    """

}