#!/usr/bin/env nextflow

/*===============
Collapse BAM file
===============*/

process BAM_COLLAPSE {

    cpus 5

    input:
    path reads

    output:
    path("*Collapsed.bam"), emit: collapsedbam
    path("*tRNAs-almost-mapped.count"), emit: tRNA_almost_mapped_count
    //path("*.log"), emit: collapse_logs

    script:
    """
    BAMcollapse.sh ${reads} ${task.cpus} > BAMcollapse.log
    """
}