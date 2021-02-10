#!/usr/bin/env nextflow

/*==========================
Calculate cleavage score
==========================*/


process COMBINED_SCORE {

    input:
    path summary
    path layout

    output:
    path("*.pdf"), emit: pdfs
    path("*absolute-score-results.tsv"), emit: abs_res_tsv
    path("*relative-score-results.tsv"), emit: rel_res_tsv
    path("*.txt"), emit: top_features
    path("*score*"), emit: all_txt_output

    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    CombinedScore.R \\
        \${condition1}_vs_\${condition2}_summarised_all-results.tsv \\
        Combined_Score
    cat Combined_Score_relative-score-results.tsv \\
        | grep -v ^combined \\
        | awk '{print \$1}' \\
        > High-combined-score-features_feature-names.txt
    """
}