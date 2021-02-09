#!/usr/bin/env nextflow

/*==========================
Calculate cleavage score
==========================*/


process FISHERS_METHOD {

    input:
    path depthfiles
    path layout
    path gtf

    output:
    path("*tsRNAs.tsv"), emit: trnas
    path("*ncRNAs.tsv"), emit: ncrnas
    path("*feature-names.txt"), emit: top_features
    path("*pvalues.tsv"), emit: pvalues
    path("*.pdf"), emit: pdfs

    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    ### Generating p-values for all features based on combination of t.tests using Fisher's method
    Pvalue_generator.R \\
        Everything_\${condition1}.depth \\
        Everything_\${condition2}.depth \\
        $gtf \\
        Feature-P-values
    sort -t\$'\t' -k2,2g \\
        Feature-P-values_FisherMethod_pvalues.tsv \\
        | grep -v Fishers.method.pvalue \\
        | awk -F '\t' '(\$2 + 0) < 0.05' \\
        > FishersMethod_Features.tsv
    grep -v ^ENS FishersMethod_Features.tsv > FishersMethod_tsRNAs.tsv
    grep ^ENS FishersMethod_Features.tsv > FishersMethod_ncRNAs.tsv
    awk -F' ' '{print \$1}' FishersMethod_Features.tsv > High-fishersmethod-score_feature-names.txt
    """
}