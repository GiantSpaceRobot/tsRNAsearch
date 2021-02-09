#!/usr/bin/env nextflow

/*==========================
Calculate cleavage score
==========================*/


process RESULTS_TABLE {

    input:
    path fishers_method
    path distribution_results
    path cleavage_results
    path deseq2_results
    path slope_results
    path layout

    output:
    path("*.tsv"), emit: tsv


    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    ResultsTable.R \\
        Feature-P-values_FisherMethod_pvalues.tsv \\
        Everything_tRNAs_cond1-vs-cond2_high-distribution-score.tsv \\
        Everything_ncRNAs_cond1-vs-cond2_high-distribution-score.tsv \\
        \${condition1}_vs_\${condition2}_tRNAs_high-cleavage-score.tsv \\
        \${condition1}_vs_\${condition2}_ncRNAs_high-cleavage-score.tsv \\
        \${condition1}_vs_\${condition2}_DESeq2-output.csv \\
        \${condition1}_vs_\${condition2}_tRNAs_high-slope-score.tsv \\
        \${condition1}_vs_\${condition2}_ncRNAs_high-slope-score.tsv \\
        \${condition1}_vs_\${condition2}
    touch my.tsv 
    """
}