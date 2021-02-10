#!/usr/bin/env nextflow

/*==========================
==========================*/

process BARPLOTS {

    input:
    path raw_counts
    path layout
    path gtf

    output:
    path("*.pdf"), emit: pdf

    script:
    """
    Barplots.R \\
        All-Features_raw-counts.tsv \\
        ${gtf} \\
        ${layout}
    """
}