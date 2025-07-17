process DATA_QC {
    
    cpus 8
    memory '50 GB'
    
    container 'quay.io/jaxcompsci/r-qtl2-deseq-biomart-tidy:v4'

    input:
    tuple val(id), path(covar_file), path(cross_file), path(genoprobs_file), path(alleleprobs_file), path(kinship_file), path(pheno_file), path(covar_info_file)

    output:
    tuple val(id), path("pr.rds"), path("apr.rds"), path("cross.rds"), path("kinship.rds"), path("covar.csv"), emit: probs_files
    path("*_pheno.csv"), emit: pheno_files
    path("*_covar_info.csv"), emit: covar_info_files

    script:

    """
    Rscript ${projectDir}/bin/qtl/data_qc.R ${covar_file} \
            ${cross_file} \
            ${genoprobs_file} \
            ${alleleprobs_file} \
            ${kinship_file} \
            ${pheno_file} \
            ${covar_info_file}
    """

}