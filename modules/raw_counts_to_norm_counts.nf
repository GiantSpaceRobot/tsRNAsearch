#!/usr/bin/env nextflow

/*==================
Normalise raw counts
==================*/

//params.output_dir = 'results'

process RAW_COUNTS_TO_NORM_COUNTS {

    tag "$countfile.simpleName"
    //publishDir params.output_dir

    input:
    path countfile
    file trna_gtf
    file ncrna_gtf
    
    output:
    path("*_tpm.count"), emit: tpm_count
    path("*_raw-tpm.tsv"), emit: tsv

    script:
    """
    Count-to-TPM.sh \\
        ${countfile} \\
        ${countfile.simpleName} \\
        ${trna_gtf} \\
        ${ncrna_gtf}
    """
}