#!/usr/bin/env nextflow

/*============================
============================*/

process CREATE_EMPTY_PDFS {

    input:
    path layout

    output:
    path("*.pdf"), emit: pdfs

    script:
    """

    """
}