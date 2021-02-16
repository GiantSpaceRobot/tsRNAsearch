#!/usr/bin/env nextflow

/*==========================
==========================*/

process PLOT_TRNAS {

    input:
    path depthfiles
    path top_features
    path layout

    output:
    path("*.pdf"), emit: pdf

    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    Depthfile_Plotter_tRNA-DEGs.R \\
        Everything_\${condition1}_concatenated_mean_stdev_tsRNA.depth \\
        Everything_\${condition2}_concatenated_mean_stdev_tsRNA.depth \\
        High-combined-score-features_feature-names.txt \\
        \${condition1}_vs_\${condition2}_Features_All-tsRNAs.pdf \\
        0 \\
        "yes" &
    Depthfile_Plotter_tRNA-DEGs.R \\
        Everything_\${condition1}_concatenated_mean_stdev_tsRNA.depth \\
        Everything_\${condition2}_concatenated_mean_stdev_tsRNA.depth \\
        High-combined-score-features_feature-names.txt \\
        \${condition1}_vs_\${condition2}_Features_Top-tsRNAs.pdf \\
        0 \\
        "no" &
    wait
    mv \${condition1}_vs_\${condition2}_Features_All-tsRNAs.pdf TranscriptionPlots_All-tsRNAs.pdf
    mv \${condition1}_vs_\${condition2}_Features_Top-tsRNAs.pdf TranscriptionPlots_Top-tsRNAs.pdf
    """
}