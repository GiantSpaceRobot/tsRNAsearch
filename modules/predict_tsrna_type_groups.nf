#!/usr/bin/env nextflow

/*=============================================================
 
=============================================================*/

process PREDICT_TSRNA_TYPE_GROUPS {

    //tag "$tsrna_depthfile.simpleName"

    input:
    path depth_means
    path layout

    output:
    path("*_HTML.txt"), emit: html
    path("*.tsv"), emit: tsv

    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    # Condition 1
    cat Combined_\${condition1}_tRNAs-almost-mapped_RPM_depth.mean \\
        tsRNA_\${condition1}_concatenated_depth.mean \\
        > \${condition1}_tsRNAs.depth.mean
    tsRNA-type-classification.R \\
        \${condition1}_tsRNAs.depth.mean \\
        \${condition1}_predicted-tsRNA-types.txt
    awk '{print \$1"\t"\$2}' \${condition1}_predicted-tsRNA-types.txt \\
        > \${condition1}_predicted-tsRNA-types_clean.txt
    # Condition 2
    cat Combined_\${condition2}_tRNAs-almost-mapped_RPM_depth.mean \\
        tsRNA_\${condition2}_concatenated_depth.mean \\
        > \${condition2}_tsRNAs.depth.mean
    tsRNA-type-classification.R \\
        \${condition2}_tsRNAs.depth.mean \\
        \${condition2}_predicted-tsRNA-types.txt
    awk '{print \$1"\t"\$2}' \${condition2}_predicted-tsRNA-types.txt \\
        > \${condition2}_predicted-tsRNA-types_clean.txt
    # Combined
    paste \${condition1}_predicted-tsRNA-types_clean.txt \\
        \${condition2}_predicted-tsRNA-types_clean.txt \\
        | grep -v ^feature \\
        > Predicted-tsRNA-types.txt
    echo -e "feature\t\${condition1}\t\${condition2}" > Predicted-tsRNA-types.tsv
    awk -F '\t' '{print \$1"\t"\$2"\t"\$4}' Predicted-tsRNA-types.txt \\
        >> Predicted-tsRNA-types.tsv
    sed -e 's/^/<br \\/>/' Predicted-tsRNA-types.tsv \
        > Predicted-tsRNA-types_HTML.txt
    """
}