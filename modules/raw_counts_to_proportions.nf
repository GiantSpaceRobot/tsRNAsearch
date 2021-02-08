#!/usr/bin/env nextflow

/*===============================

===============================*/

//params.output_dir = 'results'

process RAW_COUNTS_TO_PROPORTIONS {

    tag "$countfile.simpleName"
    publishDir params.output_dir

    input:
    path countfile
    path multimapper_countfile
    file mapped_read_counts
    //file gtf

    output:
    path("*.tsv"), emit: tsv

    script:
    """
    cat ${countfile} \\
        ${multimapper_countfile} \\
        > ${countfile.simpleName}_all-tRNAs.count
    CountCollapse.sh \\
        ${countfile.simpleName}_all-tRNAs.count \\
        ${mapped_read_counts} \\
        ${countfile.simpleName}_tRNA-mapping-information.tsv
    for i in *trimmed_accepted*; do 
        mv \$i \$(echo \$i | sed 's/_trimmed_accepted_hits_tRNAs//g')
    done
    """
}