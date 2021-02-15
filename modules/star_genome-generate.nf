#!/usr/bin/env nextflow

/*======================
Run STAR genome generate
======================*/


process MAKE_STAR_DB {

    label 'big_mem' // See nextflow.config

    input:
    path fasta

    output:
    path "star_db", emit: star_index 

    script:
    def memory  = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    echo $memory
    mkdir -p star_db
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_db/ \\
        --genomeFastaFiles $fasta \\
        --runThreadN $task.cpus \\
        --genomeSAsparseD 3 \\
        --genomeSAindexNbases 12 \\
        --genomeChrBinNbits 14 \\
        $memory
    STAR --version | sed -e "s/STAR_//g" > STAR.version.txt
    """

}
