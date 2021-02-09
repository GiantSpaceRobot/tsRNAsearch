#!/usr/bin/env nextflow

/*==========================
Calculate cleavage score
==========================*/


process CLEAVAGE_SCORE {

    input:
    path depth_means
    path layout
    path gtf

    output:
    path("*feature-names.txt"), emit: top_features
    path("*high-cleavage-score.tsv"), emit: combined_features
    path("*tRNAs*cleavage-score.tsv"), emit: tsRNA_features
    path("*ncRNAs*cleavage-score.tsv"), emit: ncRNA_features
    path("*.pdf"), emit: pdfs

    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    CleavageScore.R \\
        sorted_Everything_ncRNAs_\${condition1}_depth.mean \\
        sorted_Everything_ncRNAs_\${condition2}_depth.mean \\
        \${condition1}_vs_\${condition2}_ncRNAs \\
        $gtf &
    CleavageScore.R \\
        sorted_Everything_tRNAs_\${condition1}_depth.mean \\
        sorted_Everything_tRNAs_\${condition2}_depth.mean \\
        \${condition1}_vs_\${condition2}_tRNAs \\
        $gtf &
    wait
    cat *high-cleavage-score.tsv | grep -v ^feat | awk '{print \$1}' \\
        > High-cleavage-score_feature-names.txt
    """
}