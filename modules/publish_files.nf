#!/usr/bin/env nextflow

/*===============
Generate results PDF for groups
===============*/

process PUBLISH_FILES {

    input:
    path output_dir
    path deseq2
    path stacked_barplots
    path barplots
    path venn_diagram
    path combined_score
    path trna_plots
    path ncrna_plots
    path layout
    path deseq2_csvs
    path venn_diagram_tsvs
    path combined_score_results
 
    //output:
    //path("*.pdf"), emit: pdf
    //cp *.csv ${output_dir}/\${condition1}_vs_\${condition2}/Data/DESeq2/


    script:
    """
    condition1=\$( awk -F ',' '{print \$2}' ${layout} | uniq | head -n 1 ) # Get the first condition name
    condition2=\$( awk -F ',' '{print \$2}' ${layout} | uniq | tail -n 1 ) # Get the second condition name
    mkdir -p ${output_dir}/\${condition1}_vs_\${condition2}
    mkdir -p ${output_dir}/\${condition1}_vs_\${condition2}/Plots
    mkdir -p ${output_dir}/\${condition1}_vs_\${condition2}/Data
    mkdir -p ${output_dir}/\${condition1}_vs_\${condition2}/Data/VennDiagram
    mkdir -p ${output_dir}/\${condition1}_vs_\${condition2}/Data/Distribution
    mkdir -p ${output_dir}/\${condition1}_vs_\${condition2}/Data/Slope
    mkdir -p ${output_dir}/\${condition1}_vs_\${condition2}/Data/Cleavage
    mkdir -p ${output_dir}/\${condition1}_vs_\${condition2}/Data/FishersMethod
    mkdir -p ${output_dir}/\${condition1}_vs_\${condition2}/Data/DESeq2

    cp *intersect*.tsv ${output_dir}/\${condition1}_vs_\${condition2}/Data/VennDiagram/
    cp Combined*.tsv ${output_dir}/\${condition1}_vs_\${condition2}/Data/
    cp *.pdf ${output_dir}/\${condition1}_vs_\${condition2}/Plots/

    """
}