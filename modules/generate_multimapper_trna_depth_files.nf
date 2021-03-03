#!/usr/bin/env nextflow

/*========================================
Generate depth files for multimapper reads
========================================*/

params.species = "human"  // Default species. This updates using main.nf parameters

process GENERATE_MULTIMAPPER_TRNA_DEPTH_FILES {

    tag "$txtfile.simpleName"

    input:
    file txtfile
    file trna_gtf
    file ncrna_gtf

    output:
    path '*sorted.depth', emit: depth_files

    script:
    """
    Multimappers-to-Depthfile.py \\
        ${txtfile} \\
        "$projectDir"/additional-files/"$params.species"_tRNA-lengths.txt \\
        "${txtfile.simpleName}".depth
    NormaliseDepthFile_TPM.sh \\
        "$txtfile.simpleName".depth \\
        "$txtfile.simpleName"_normalised.depth \\
        ${trna_gtf} \\
        ${ncrna_gtf}
    sort -k1,1 -k2,2n "$txtfile.simpleName"_normalised.depth \\
        > "$txtfile.simpleName"_sorted.depth
    """

}