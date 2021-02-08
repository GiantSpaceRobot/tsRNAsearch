#!/usr/bin/env nextflow

/*===============
Collapse BAM file
===============*/

//task.cpus = 'results'

process BAM_COLLAPSE {

    cpus 20 // Max cpus usable by all processes or 20 cpus are used in each process regardless or other running processes? I think the former but not sure

    input:
    path reads

    output:
    path("*Collapsed.bam"), emit: collapsedbam
    path("*tRNAs-almost-mapped.count"), emit: tRNA_almost_mapped_count
    path("*tRNAs-almost-mapped.txt"), emit: tRNA_almost_mapped_txt

    //path("*.log"), emit: collapse_logs

    script:
    """
    BAMcollapse.sh ${reads} ${task.cpus} > BAMcollapse.log
    """
}