#!/usr/bin/env nextflow

/*============================
Split layout file by condition
============================*/

process SPLIT_LAYOUT {

    tag "$layout.simpleName"

    input:
    path layout

    output:
    val condition1, emit: cond1
    val condition2, emit: cond2

    script:
    """
    condition1=$( awk -F ',' '{print \$2}' "\$layout" | uniq | head -n 1 ) # Get the first condition name
    condition2=$( awk -F ',' '{print \$2}' "\$layout" | uniq | tail -n 1 ) # Get the second condition name
    """
}