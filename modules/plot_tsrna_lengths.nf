#!/usr/bin/env nextflow

/*==========================
==========================*/


process PLOT_TSRNA_LENGTHS {

    input:
    path tsrna_lengths
    path layout

    output:
    path("*tRNA-read-alignment-lengths.txt"), emit: txts
    path("*.pdf"), emit: pdfs

    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    grep \$condition1 ${layout} > condition1.csv
    grep \$condition2 ${layout} > condition2.csv
    # Condition 1
    mkdir -p \$condition1
    cat condition1.csv | while read i
    do
        echo \$i
        sample_name=\$(echo \$i | awk -F '.' '{print \$1}')
        echo \$sample_name
        cp \${sample_name}_tRNA-alignment-length.txt \$condition1/
    done
    # Condition 2
    mkdir -p \$condition2
    cat condition2.csv | while read i
    do
        echo \$i
        sample_name=\$(echo \$i | awk -F '.' '{print \$1}')
        echo \$sample_name
        cp \${sample_name}_tRNA-alignment-length.txt \$condition2/
    done
    tRNA_Alignment_Length_multi-replicates.R \\
        \$condition1/ \\
        \${condition1}
    tRNA_Alignment_Length_multi-replicates.R \\
        \$condition2/ \\
        \${condition2}
    """
}