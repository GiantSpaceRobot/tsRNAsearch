#!/usr/bin/env nextflow

/*==========================
==========================*/

process STACKED_BARPLOTS {

    input:
    path tsvs

    output:
    path("*.pdf"), emit: pdf

    script:
    """
    StackedBarplots.R \\
        ./ \\
        Stacked-Barplot   
    """
}