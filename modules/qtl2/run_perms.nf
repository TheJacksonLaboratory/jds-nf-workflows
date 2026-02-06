process RUN_PERMS {

    tag "$phenotype"
    
    cpus 32
    memory '200 GB'
    time 24.hour
    
    container 'quay.io/jaxcompsci/r-qtl2-deseq-biomart-tidy:v4'

    publishDir "${params.pubdir}/${id}/${phenotype}", pattern:"*_scan1perms.rds", mode:'copy'

    input:
    tuple val(id), path(genoprobs_file), path(alleleprobs_file), path(kinship_file), path(covar_file), val(phenotype), path(pheno_file), path(covar_info_file), path(map_file)

    output:
    tuple val(id), val(phenotype), path('*_scan1perms.rds'), emit: perm_files

    script:

    """
    Rscript ${projectDir}/bin/qtl/run_perms.R ${covar_file} \
            ${map_file} \
            ${genoprobs_file} \
            ${alleleprobs_file} \
            ${kinship_file} \
            ${pheno_file} \
            ${covar_info_file} \
            ${params.n_perms} \
            ${task.cpus}
    """

}