#!/usr/bin/env nextflow

/*==========================
==========================*/

process PLOT_ALL_NCRNAS {

    input:
    path depthfiles
    path layout
    path gtf

    output:
    path("*.pdf"), emit: pdfs

    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    echo Hello > non-empty_file.txt
    Depthfile_Plotter_DEGs.R \\
        \${condition1}_concatenated_mean_stdev_ncRNA.depth \\
        \${condition2}_concatenated_mean_stdev_ncRNA.depth \\
        non-empty_file.txt \\
        \${condition1}_vs_\${condition2}_Features_All-ncRNAs.pdf \\
        0 \\
        "yes" \\
        $gtf
    mv \${condition1}_vs_\${condition2}_Features_All-ncRNAs.pdf  TranscriptionPlots_All-ncRNAs.pdf
    """
}