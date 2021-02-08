#!/usr/bin/env nextflow

/*==================
Normalise raw counts
==================*/

//params.output_dir = 'results'

process RAW_COUNTS_TO_NORM_COUNTS {

    tag "$countfile.simpleName"
    //publishDir params.output_dir

    input:
    path countfile
    file mapped_read_counts
    //file gtf

    output:
    path("*_rpm.count"), emit: rpm_count
    path("*_raw-rpm.tsv"), emit: tsv

    script:
    """
    Count-to-RPM.sh \\
        ${countfile} \\
        ${mapped_read_counts} \\
        ${countfile.simpleName}
    """
}