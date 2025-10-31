process HARVEST_QTL {

    tag "$id"
    
    cpus 8
    memory '50 GB'
    
    container 'quay.io/jaxcompsci/r-qtl2-deseq-biomart-tidy:v4'

    publishDir "${params.pubdir}/${id}/qtl_scans", pattern:"*_scan1_thresh.png", mode:'copy'

    input:
    tuple val(id), path(perm_files), val(phenos), path(scan1_files), path(map_files)

    output:
    tuple val(id), path("peaks.csv"), emit: qtl_table
    tuple val(id), path(map_files), path(scan1_files), emit: qtl_files
    tuple val(id), path("*_scan1_thresh.png"), emit: scan1_thresh

    script:

    """
    Rscript ${projectDir}/bin/qtl/harvest_qtl.R
    """
}