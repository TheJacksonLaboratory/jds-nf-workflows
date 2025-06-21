process MAP_QTL {

    tag "$phenotype"
    
    cpus 8
    memory '50 GB'
    
    container 'quay.io/jaxcompsci/r-qtl2-deseq-biomart-tidy:v4'

    publishDir "${params.pubdir}/${id}", pattern:"*_scan1.png", mode:'copy'
    publishDir "${params.pubdir}/${id}", pattern:"*_scan1out.rds", mode:'copy'

    input:
    tuple val(id), path(genoprobs_file), path(alleleprobs_file), path(cross_file), path(kinship_file), path(covar_file), path(pheno_file), val(phenotype)

    output:
    tuple val(id), path(genoprobs_file), path(alleleprobs_file), path("*_cross.rds"), path(kinship_file), path(covar_file), path(pheno_file), val(phenotype), emit: probs_files
    tuple val(id), path("*_scan1out.rds"), path("*_cross.rds"), emit: scan1_files
    tuple val(id), path("*_scan1.png"), emit: scan1_plots

    script:

    """
    Rscript ${projectDir}/bin/qtl/map_qtl.R ${covar_file} \
            ${cross_file} \
            ${genoprobs_file} \
            ${alleleprobs_file} \
            ${kinship_file} \
            ${pheno_file}

    mv ${cross_file} ${phenotype}_cross.rds
    """

}