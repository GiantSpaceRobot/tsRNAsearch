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
        ${depthfile.simpleName}_tRNAs_Transcription-plots.pdf \\
        0
    Single-replicate-analysis.R \\
        ${depthfile} \\
        ${depthfile.simpleName}
    # If there are less than 21 lines, get max lines
    if [[ \$(wc -l <${depthfile.simpleName}_high-combined-score.tsv) -ge 21 ]];
    then
        head_lines=21
        tail_lines=20
    else
        # Head = max file lines. Tail = head - 1
        head_lines=\$(wc -l <${depthfile.simpleName}_high-combined-score.tsv)
        tail_lines="\$((\$head_lines-1))"
    fi
    head -n \$head_lines ${depthfile.simpleName}_high-combined-score.tsv | tail -n \$tail_lines | while read line
    do
        top_tRNA=\$(echo \$line | awk '{print \$1}')
        grep -w "\$top_tRNA" ${depthfile} >> ${depthfile.simpleName}_top-20-tRNAs.depth
    done
    Depthfile_Plotter.R \\
        ${depthfile.simpleName}_top-20-tRNAs.depth \\
        ${depthfile.simpleName}_top-20-tRNAs_Transcription-Plots.pdf \\
        0
    for i in *trimmed_accepted*; do 
        mv \$i \$(echo \$i | sed 's/_trimmed_accepted_hits_tRNAs_sorted//g')
    done
    """
}