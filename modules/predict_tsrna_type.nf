#!/usr/bin/env nextflow

/*=============================================================
 
=============================================================*/

process PREDICT_TSRNA_TYPE {

    //tag "$tsrna_depthfile.simpleName"

    input:
    path tsrna_depthfile
    path multimapper_tsrna_depthfile

    output:
    path("*_HTML.txt"), emit: html
    path("*_RmdHTML.txt"), emit: rmkdwn_html

    shell:
    '''
    echo !{tsrna_depthfile} | tr " " "\n" | sort > tsRNA_Depthfiles.txt
    Predict-tsRNA-Types.sh tsRNA_Depthfiles.txt > Predict-tsRNA-Types.log
    '''
}