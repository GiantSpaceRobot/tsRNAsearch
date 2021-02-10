#!/usr/bin/env nextflow

/*==========================
==========================*/

process GENERATE_COUNT_DATAFRAME {

    input:
    path countfiles

    output:
    path("*.tsv"), emit: count_df

    script:
    """
    ### RPM count files
    # Get feature names from one of the count files
    cat `ls *collapsed.count | head -n 1` \\
        | awk '{print \$1}' \\
        > FCount.rpm.all-features
    echo "Feature" > Count_Header_File.txt
    for f in *rpm_collapsed.count; do
        echo \$f
        file_base=\$(basename \$f)
        filename="\$( cut -d '.' -f 1 <<< "\$file_base" )"
        my_header=\$(cat Count_Header_File.txt)
        new_header=\$(echo -e \$my_header"\t"\$filename)
        echo \$new_header > Count_Header_File.temp.txt
        mv Count_Header_File.temp.txt Count_Header_File.txt
        awk '{print \$2}' \$f | paste FCount.rpm.all-features - \\
            >> FCount.rpm.temp
        mv FCount.rpm.temp FCount.rpm.all-features
    done
    cat Count_Header_File.txt FCount.rpm.all-features \\
        > All-Features_normalised-counts.tsv
    
    ### Raw count files
    # Get feature names from one of the count files
    cat `ls *counts_collapsed.count | head -n 1` \\
        | awk '{print \$1}' \\
        > FCount.raw.all-features
    echo "Feature" > Count_Header_File.txt
    for f in *counts_collapsed.count; do
        echo \$f
        file_base=\$(basename \$f)
        filename="\$( cut -d '.' -f 1 <<< "\$file_base" )"
        my_header=\$(cat Count_Header_File.txt)
        new_header=\$(echo -e \$my_header"\t"\$filename)
        echo \$new_header > Count_Header_File.temp.txt
        mv Count_Header_File.temp.txt Count_Header_File.txt
        awk '{print \$2}' \$f | paste FCount.raw.all-features - \\
            >> FCount.raw.temp
        mv FCount.raw.temp FCount.raw.all-features
    done
    cat Count_Header_File.txt FCount.raw.all-features \
        > All-Features_raw-counts.tsv
    """
}