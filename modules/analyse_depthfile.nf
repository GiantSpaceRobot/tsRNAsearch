#!/usr/bin/env nextflow

/*=============================================================
Geenrate mean and standard deviation TSV files from depth files
=============================================================*/

process ANALYSE_DEPTHFILE {

    tag "$depthfile.simpleName"

    input:
    path depthfile

    output:
    path("*_mean-std_sorted.tsv"), emit: stats
    path("*_depth_stats.tsv"), emit: depth_inf

    script:
    """
    Analyse-depthfile.py \\
        $depthfile \\
        "$depthfile.simpleName".tsv
    Sort-depthfile-tsv.sh \\
        "$depthfile.simpleName".tsv \\
        "$depthfile.simpleName"_mean-std_sorted.tsv
    Generate-depthfile-stats.sh \\
        "$depthfile.simpleName"_mean-std_sorted.tsv \\
    """
}