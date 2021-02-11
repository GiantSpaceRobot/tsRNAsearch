#!/usr/bin/env nextflow

/*===============
Move all results files for each sample into a directory of its own
===============*/

process ORGANISE_RESULTS_GROUPS {

    input:
    path output_dir
    path mapped_read_counts
    path pdfs

    output:
    path("*finished"), emit: signal


    script:
    """
    cat ${mapped_read_counts} | grep -v ^Sample | while read line
    do
        my_sample=\$(echo \$line | awk -F '_trimmed' '{print \$1}')
        echo \$my_sample
        mkdir -p ${output_dir}/temp_dir
        mv ${output_dir}/*\$my_sample* ${output_dir}/temp_dir/
        mkdir -p ${output_dir}/\$my_sample
        mv ${output_dir}/temp_dir/* ${output_dir}/\$my_sample/
        rmdir ${output_dir}/temp_dir
    done
    touch finished
    """
}