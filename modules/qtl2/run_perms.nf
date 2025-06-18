process RUN_PERMS {

    tag "$phenotype"
    
    cpus 8
    memory '50 GB'
    
    container 'quay.io/jaxcompsci/r-qtl2-deseq-biomart-tidy:v4'

    input:
    tuple val(id), path(genoprobs_file), path(alleleprobs_file), path(cross_file), path(kinship_file), path(covar_file), path(pheno_file), val(phenotype)

    output:
    tuple val(id), path(genoprobs_file), path(alleleprobs_file), path(cross_file), path(kinship_file), path(covar_file), path(pheno_file), val(phenotype), emit: probs_files
    tuple val(id), val(phenotype), path("*_scan1perms.rds"), emit: perm_files

    script:

    """
    Rscript ${projectDir}/bin/qtl/run_perms.R ${covar_file} \
            ${cross_file} \
            ${genoprobs_file} \
            ${alleleprobs_file} \
            ${kinship_file} \
            ${pheno_file} \
            ${params.n_perms}

    """

}