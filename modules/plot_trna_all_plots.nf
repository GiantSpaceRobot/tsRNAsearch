#!/usr/bin/env nextflow

/*===============================
Analyse tRNAs from single samples
===============================*/

params.output_dir = 'results'

process PLOT_TRNA_ALL_PLOTS {

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
    Depthfile_Plotter.R \\
        ${depthfile} \\
        ${depthfile.simpleName}_tRNAs_Coverage-plots.pdf \\
        0 \\
        ${gtf}
    Single-replicate-analysis.R \\
        ${depthfile} \\
        ${depthfile.simpleName}_tRNAs \\
        ${gtf}
    for i in *trimmed_accepted*; do 
        mv \$i \$(echo \$i | sed 's/_trimmed_accepted_hits_tRNAs_sorted//g')
    done
    """
}