#!/usr/bin/env nextflow

/*========
Run MultiQC
========*/

process MULTIQC {

    input:
    file ('fastqc/*')

    output:
    file "*multiqc_report.html"
    file "*_data"

    script:
    """
    multiqc . -m fastqc
    """
}