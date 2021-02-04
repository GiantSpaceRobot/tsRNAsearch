#!/usr/bin/env nextflow

/*=================================
Generate depth files from BAM files
=================================*/

params.species = "human"  // Default species. This updates using main.nf parameters

process GENERATE_TRNA_DEPTH_FILES {

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
    Remove-introns.R \\
        "$bamfile.simpleName"_normalised.depth \\
        "$projectDir"/additional-files/"$params.species"_tRNA-introns-for-removal.tsv \\
        "$bamfile.simpleName"_introns-removed.depth
    Bedgraph_collapse-tRNAs.py \\
        "$bamfile.simpleName"_introns-removed.depth \\
        "$bamfile.simpleName"_collapsed.depth
    sort -k1,1 -k2,2n "$bamfile.simpleName"_collapsed.depth \\
        > "$bamfile.simpleName"_sorted.depth
    """

}