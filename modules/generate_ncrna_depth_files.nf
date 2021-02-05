#!/usr/bin/env nextflow

/*=================================
Generate depth files from BAM files
=================================*/

params.species = "human"  // Default species. This updates using main.nf parameters

process GENERATE_NCRNA_DEPTH_FILES {

    tag "$bamfile.simpleName"

    input:
    file bamfile
    file mapped_read_counts

    output:
    path '*sorted.depth', emit: depth_files

    script:
    """
    samtools depth \\
        -d 100000000 \\
        -aa $bamfile \\
        > "$bamfile.simpleName"_raw.depth
    NormaliseDepthFile.sh \\
        "$bamfile.simpleName"_raw.depth \\
        "$mapped_read_counts" \\
        "$bamfile.simpleName"_normalised.depth
    sort -k1,1 -k2,2n "$bamfile.simpleName"_normalised.depth \\
        > "$bamfile.simpleName"_sorted.depth
    """
}