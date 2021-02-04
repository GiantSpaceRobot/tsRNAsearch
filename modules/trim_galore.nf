#!/usr/bin/env nextflow

/*========
Trim reads
==========*/

process TRIM_READS {

    input:
    file reads
    val min_read_length

    output:
    path '*_trimming_report.txt', emit: trim_reports
    path '*_trimmed.fq.gz', emit: trimmed_reads

    script:
    """
    trim_galore \\
        --stringency 10 \\
        --length $min_read_length \\
        $reads
    """

}