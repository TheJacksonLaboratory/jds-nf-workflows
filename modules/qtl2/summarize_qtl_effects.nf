process SUMMARIZE_QTL_EFFECTS {

    tag "$id"
    
    cpus 8
    memory '50 GB'
    
    container 'quay.io/jaxcompsci/r-qtl2-deseq-biomart-tidy:v4'

    publishDir "${params.pubdir}/${id}", pattern:"*_qtl_effects_heatmap.png", mode:'copy'
    publishDir "${params.pubdir}/${id}", pattern:"*_all_peaks_file.csv", mode:'copy'
    publishDir "${params.pubdir}/${id}", pattern:"*_peaks_viewer_file.csv", mode:'copy'

    input:
    tuple val(id), path(peak_files)
    
    output:
    tuple val(id), path("*_qtl_effects_heatmap.png"), path("*_peaks_viewer_file.csv"), path("*_all_peaks_file.csv"), emit: qtl_effects_summary_files
    

    script:

    """
    Rscript ${projectDir}/bin/qtl/summarize_qtl_effects.R ${id} ${params.primary_chrom_bed}
    """

}