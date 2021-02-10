#!/usr/bin/env nextflow

/*=================
Collapse counts
=================*/


process COUNTS_TO_COLLAPSED_COUNTS {

    tag "$countfile.simpleName"

    input:
    path countfile

    output:
    path("*collapsed.count"), emit: collapsed_count

    script:
    """
    CollapseCountfile.py \\
        ${countfile} \\
        ${countfile.simpleName}_collapsed.count
    """
}