#!/usr/bin/env nextflow

/*================================
Analyse ncRNAs from single samples
================================*/

params.output_dir = 'results'

process PLOT_NCRNA_ALL_PLOTS {

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
    #for i in \$(seq 1 5); do touch ${depthfile.simpleName}_\$i.pdf; done
    Depthfile_Plotter.R \\
        ${depthfile} \\
        ${depthfile.simpleName}_ncRNAs_Coverage-plots.pdf \\
        0 \\
        ${gtf}
    Single-replicate-analysis.R \\
        ${depthfile} \\
        ${depthfile.simpleName}_ncRNAs \\
        ${gtf}
    for i in *trimmed_accepted*; do 
        mv \$i \$(echo \$i | sed 's/_trimmed_accepted_hits_ncRNAs_sorted//g')
    done
    """
}