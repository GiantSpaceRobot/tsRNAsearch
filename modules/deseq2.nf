#!/usr/bin/env nextflow

/*=================================================
Carry out differential gene expression using DESeq2
=================================================*/

process DESEQ2 {

    input:
    path countfiles
    path layout
    path gtf

    output:
    path("*.pdf"), emit: pdfs
    path("*.csv"), emit: csvs

    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    echo Conditions are \$condition1 and \$condition2
    DESeq2_tsRNAsearch.R ./ "\${condition1}_vs_\${condition2}" ${gtf} ${layout}
    """
}