#!/usr/bin/env nextflow

/*=================================
Generate depth files from BAM files
=================================*/

params.species = "human"  // Default species. This updates using main.nf parameters

process GENERATE_DEPTH_FILES {

    input:
    file bamfile

    output:
    path '*.depth', emit: depth_files

    script:
    """
    samtools depth \\
        -d 100000000 \\
        -aa $bamfile \\
        > "$bamfile.simpleName"_raw.depth
    Rscript bin/Remove-introns.R \\
        "$bamfile.simpleName"_raw.depth \\
        "$projectDir/additional-files/${species}_tRNA-introns-for-removal.tsv" \\
        "$bamfile.simpleName"_introns-removed.depth
    

    """

}