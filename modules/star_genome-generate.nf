#!/usr/bin/env nextflow

/*======================
Run STAR genome generate
======================*/


process MAKE_STAR_DB {

    input:
    path fasta

    output:
    path "star_db", emit: star_index 

    """
    mkdir -p star_db
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_db/ \\
        --genomeFastaFiles $fasta \\
        --genomeSAindexNbases 8 \\
        --runThreadN $task.cpus
    STAR --version | sed -e "s/STAR_//g" > STAR.version.txt
    """

}
