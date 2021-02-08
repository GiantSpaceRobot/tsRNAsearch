#!/usr/bin/env nextflow

/*===============
Generate results PDF
===============*/

params.output_dir = 'results'

process GENERATE_RESULTS_PDF {

    publishDir params.output_dir

    input:
    path alignment_pdf
    path analysis_pdfs
    path mapped_read_counts

    output:
    path("*Results-summary.pdf"), emit: pdf

    script:
    """
    cat ${mapped_read_counts} | grep -v ^Sample | while read line
    do
        my_sample=\$(echo \$line | awk -F '_trimmed' '{print \$1}')
        echo Creating results PDF for \$my_sample
        mkdir \$my_sample
        cp *.pdf \$my_sample/
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dAutoRotatePages=/None \\
            -sOutputFile=\${my_sample}_Results-summary.pdf \\
            "$projectDir/additional-files/tsRNAsearch_page.pdf" \\
            \${my_sample}_tRNA-alignment-length.pdf \\
            \${my_sample}_high-distribution-score.pdf \\
            \${my_sample}_high-cleavage-score.pdf \\
            \${my_sample}_high-slope-score.pdf \\
            \${my_sample}_high-combined-score.pdf \\
            \${my_sample}_top-20-tRNAs_Transcription-Plots.pdf
    done
    """
}