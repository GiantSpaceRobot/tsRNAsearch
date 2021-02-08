#!/usr/bin/env nextflow

/*=====================================================================
Expand multimapper count files by adding tRNA families with counts of 0
=====================================================================*/

//params.species = 'human'

process ADD_EMPTY_COUNTS {

    tag "$multi_count.simpleName"

    input:
    path multi_count
    path empty_count

    output:
    path("*filled.count"), emit: filled_multi_counts

    script:
    """
    FillEmptyCounts.py \\
        ${multi_count} \\
        ${empty_count} \\
        ${multi_count.simpleName}_filled_unsorted.count
    # Sort and remove empty lines
    sort -k1,1 ${multi_count.simpleName}_filled_unsorted.count > ${multi_count.simpleName}_filled.count
    """
}