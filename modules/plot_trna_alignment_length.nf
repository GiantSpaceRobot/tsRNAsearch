#!/usr/bin/env nextflow

/*=========================
Plot tRNA alignment lengths
=========================*/

params.output_dir = 'results'

process PLOT_TRNA_ALIGNMENT_LENGTH {

    tag "$bamfile.simpleName"
    publishDir params.output_dir

    input:
    path bamfile

    output:
    path("*tRNA-alignment-length.pdf"), emit: pdf
    path("*tRNA-alignment-length.txt"), emit: txt

    script:
    """
    samtools view -h -o ${bamfile.simpleName}.sam ${bamfile}
    grep -v ^@ ${bamfile.simpleName}.sam > ${bamfile.simpleName}_no-header.sam
    tRNA_Alignment_Length.R \\
        ${bamfile.simpleName}_no-header.sam \\
        ${bamfile.simpleName}_tRNA-alignment-length
    rm ${bamfile.simpleName}.sam ${bamfile.simpleName}_no-header.sam
    for i in *trimmed_accepted*; do 
        mv \$i \$(echo \$i | sed 's/_trimmed_accepted_hits_tRNAs//g')
    done
    """
}