#!/usr/bin/env nextflow

/*==========================
Generate Venn diagram
==========================*/

process VENN_DIAGRAM {

    input:
    path fishers_method
    path distribution_results
    path cleavage_results
    path deseq2_results
    path slope_results
    path layout

    output:
    path("*.pdf"), emit: pdf
    path("*intersect*.tsv"), emit: tsvs

    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    VennDiagram.R \\
        DEGs_names-only.txt \\
        High-distribution-scores_feature-names.txt \\
        High-cleavage-score_feature-names.txt \\
        High-fishersmethod-score_feature-names.txt \\
        High-slope-score_feature-names.txt \\
        \${condition1}_vs_\${condition2}
    """
}