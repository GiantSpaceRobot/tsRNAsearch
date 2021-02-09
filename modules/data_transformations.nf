#!/usr/bin/env nextflow

/*=================
Collapse raw counts
=================*/


process DATA_TRANSFORMATIONS {

    input:
    path layout
    path trna_depthfile
    path ncrna_depthfile
    path multi_trna_depthfile
    path mapped_read_counts

    output:
    path("*depth.mean"), emit: depth_means
    path("Everything*depth"), emit: depthfiles
    path("*Everything_tRNAs.stddev"), emit: trna_stddev
    path("*Everything_ncRNAs.stddev"), emit: ncrna_stddev

    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    grep -w "\$condition1" ${layout} > Condtion1.csv
    grep -w "\$condition2" ${layout} > Condtion2.csv
    GatherDepthFiles.sh Condtion1.csv \$condition1
    GatherDepthFiles.sh Condtion2.csv \$condition2
    paste \\
        Everything_\${condition1}_depth.mean \\
        Everything_\${condition2}_depth.mean \\
        > EverythingTest_\${condition1}-vs-\${condition2}.mean
    ### Calculate StdDev
    Mean-to-RelativeDifference.py \\
        EverythingTest_\${condition1}-vs-\${condition2}.mean \\
        EverythingTest_\${condition1}-vs-\${condition2}.stddev
    ### Separate into tRNAs and other ncRNAs
    grep ^ENS EverythingTest_\${condition1}-vs-\${condition2}.stddev \\
        > \${condition1}-vs-\${condition2}_Everything_ncRNAs.stddev
    grep -v ^ENS EverythingTest_\${condition1}-vs-\${condition2}.stddev \\
        > \${condition1}-vs-\${condition2}_Everything_tRNAs.stddev
    """
}