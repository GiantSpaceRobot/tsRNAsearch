#!/usr/bin/env nextflow

/*========================================
Generate depth files for multimapper reads
========================================*/

params.species = "human"  // Default species. This updates using main.nf parameters

process GENERATE_MULTIMAPPER_TRNA_DEPTH_FILES {

    tag "$txtfile.simpleName"

    input:
    file txtfile
    file mapped_read_counts

    output:
    path '*sorted.depth', emit: depth_files

    script:
    """
    Multimappers-to-Depthfile.py \\
        ${txtfile} \\
        "$projectDir"/additional-files/"$params.species"_tRNA-lengths.txt \\
        "${txtfile.simpleName}".depth
    NormaliseDepthFile.sh \\
        "$txtfile.simpleName".depth \\
        "$mapped_read_counts" \\
        "$txtfile.simpleName"_normalised.depth
    sort -k1,1 -k2,2n "$txtfile.simpleName"_normalised.depth \\
        > "$txtfile.simpleName"_sorted.depth
    """

}