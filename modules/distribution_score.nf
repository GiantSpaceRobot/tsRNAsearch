#!/usr/bin/env nextflow

/*==========================
Calculate distribution score
==========================*/


process DISTRIBUTION_SCORE {

    input:
    path ncrna_stddev
    path trna_stddev
    path gtf

    output:
    path("*feature-names.txt"), emit: top_features
    path("Everything*high-distribution-score.tsv"), emit: combined_features
    path("*tRNAs*distribution-score.tsv"), emit: tsRNA_features
    path("*ncRNAs*distribution-score.tsv"), emit: ncRNA_features
    path("*.pdf"), emit: pdfs

    script:
    """
    DistributionScore.sh $ncrna_stddev $trna_stddev $gtf
    """
}