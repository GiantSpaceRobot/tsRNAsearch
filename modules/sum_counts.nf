#!/usr/bin/env nextflow

/*============================
Sum all counts for each sample
============================*/

process SUM_COUNTS {

    input:
    file trna_counts
    file ncrna_counts
    file almost_mapped_counts

    output:
    path("*Mapped.txt"), emit: sum_counts
   
    shell:
    '''
    echo !{trna_counts} | tr " " "\n" | sort > tRNA_Counts.txt
    echo !{ncrna_counts} | tr " " "\n" | sort > ncRNA_Counts.txt
    echo !{almost_mapped_counts} | tr " " "\n" | sort > Almost_Mapped_Counts.txt
    CombineCountFiles.sh tRNA_Counts.txt > CombineFiles.log
    echo -e "Sample\tMapped_reads" > Reads_Mapped.txt
    cat *_reads-mapped.txt | sort >> Reads_Mapped.txt
    '''
}