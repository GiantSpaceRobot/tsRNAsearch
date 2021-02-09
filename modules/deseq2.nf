#!/usr/bin/env nextflow

/*=================================================
Carry out differential gene expression using DESeq2
=================================================*/

process DESEQ2 {

    input:
    path countfiles
    path layout
    path gtf

    output:
    path("*.pdf"), emit: pdfs
    path("*.csv"), emit: csvs
    path("*names-only.txt"), emit: txt


    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    echo Conditions are \$condition1 and \$condition2
    DESeq2_tsRNAsearch.R ./ "\${condition1}_vs_\${condition2}" ${gtf} ${layout}
    # Create empty results files if none were generated
    touch -a \${condition1}_vs_\${condition2}_DESeq2-output.csv
    touch -a \${condition1}_vs_\${condition2}_DESeq2-output-upregulated.csv
    touch -a \${condition1}_vs_\${condition2}_DESeq2-output-downregulated.csv
    cat *regulated.csv | grep -v ^, | sort -t',' -k7,7g | awk -F ',' '{print \$1}' | awk -F ' ' '{print \$1}' \\
        > DEGs_names-only.txt
    touch -a DEGs_names-only.txt
    """
}