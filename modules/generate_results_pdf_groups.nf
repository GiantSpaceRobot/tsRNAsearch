#!/usr/bin/env nextflow

/*===============
Generate results PDF for groups
===============*/

params.output_dir = 'results'

process GENERATE_RESULTS_PDF_GROUPS {

    publishDir params.output_dir

    input:
    path deseq2
    path stacked_barplots
    path barplots
    path venn_diagram
    path combined_score
    path trna_plots
    path ncrna_plots
    path layout

    output:
    path("*.pdf"), emit: pdf

    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    echo Creating results PDF
    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dAutoRotatePages=/None \\
        -sOutputFile=\${condition1}_vs_\${condition2}_Results-summary.pdf \\
        "$projectDir/additional-files/tsRNAsearch_page.pdf" \\
        \${condition1}_vs_\${condition2}_RPM-PCA.pdf \\
        \${condition1}_vs_\${condition2}_Distance-Matrix.pdf \\
        \${condition1}_vs_\${condition2}_BarPlot_Raw-readcounts.pdf \\
        \${condition1}_vs_\${condition2}_BarPlot_RPM-normalised.pdf \\
        Barplot_All-Samples.pdf \\
        Stacked-Barplot_tRNA-species-per-sample.pdf \\
        \${condition1}_vs_\${condition2}_VennDiagram.pdf \\
        "$projectDir/additional-files/tsRNAsearch_Combined_page.pdf" \\
        Combined_Score_tsRNAs.pdf \\
        TranscriptionPlots_Top-tsRNAs.pdf \\
        Combined_Score_ncRNAs.pdf \\
        TranscriptionPlots_Top-ncRNAs.pdf
    """
}