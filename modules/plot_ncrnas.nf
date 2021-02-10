#!/usr/bin/env nextflow

/*==========================
Plot top ncRNA features identified using combined score
==========================*/

process PLOT_NCRNAS {

    input:
    path depthfiles
    path top_features
    path layout
    path gtf

    output:
    path("*.pdf"), emit: pdfs

    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    Depthfile_Plotter_DEGs.R \\
        \${condition1}_concatenated_mean_stdev_ncRNA.depth \\
        \${condition2}_concatenated_mean_stdev_ncRNA.depth \\
        High-combined-score-features_feature-names.txt \\
        \${condition1}_vs_\${condition2}_Features_Top-ncRNAs.pdf \\
        0 \\
        "no" \\
        $gtf
    mv \${condition1}_vs_\${condition2}_Features_Top-ncRNAs.pdf  TranscriptionPlots_Top-ncRNAs.pdf
    """
}