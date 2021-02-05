#!/usr/bin/env nextflow

/*===============================
Analyse tRNAs from single samples
===============================*/

params.output_dir = 'results'

process PLOT_TRNA_SINGLE_SAMPLE_ANALYSIS {

    tag "$depthfile.simpleName"
    publishDir params.output_dir

    input:
    path depthfile
    file gtf

    output:
    path("*.pdf"), emit: pdfs
    path("*.tsv"), emit: tsvs

    script:
    """
    Single-replicate-analysis.R \\
        ${depthfile} \\
        ${depthfile.simpleName}_tRNAs \\
        ${gtf}
    for i in *trimmed_accepted*; do 
        mv \$i \$(echo \$i | sed 's/_trimmed_accepted_hits_tRNAs_sorted//g')
    done
    """
}