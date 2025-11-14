process QTL_EFFECTS {

    tag "$phenotype"
    
    cpus 8
    memory '50 GB'
    
    container 'quay.io/jaxcompsci/r-qtl2-deseq-biomart-tidy:v4'

    publishDir "${params.pubdir}/${id}/${phenotype}", pattern:"*_scan1coef.png", mode:'copy'
    publishDir "${params.pubdir}/${id}/${phenotype}", pattern:"*_scan1blup.png", mode:'copy'
    publishDir "${params.pubdir}/${id}/${phenotype}", pattern:"*_qtl_effect_files.RData", mode:'copy'
    
    input:
    tuple val(id), val(phenotype), path(genoprobs_file), path(alleleprobs_file), path(kinship_file), path(covar_file), path(pheno_file), path(covar_info_file), path(map_file), val(chrom), val(peak_pos), val(start_pos), val(end_pos), path(scan1out_file)

    output:
    tuple val(id), path("*_peaks.csv"), emit: qtl_peaks_files
    tuple val(id), path("*_qtl_effect_files.RData"), emit: qtl_effects_files
    tuple val(id), path("*_scan1coef.png"), path("*_scan1blup.png"), emit: qtl_effects_plots
    

    script:

    """
    Rscript ${projectDir}/bin/qtl/qtl_effects.R ${phenotype} \
        ${alleleprobs_file} \
        ${kinship_file} \
        ${covar_file} \
        ${pheno_file} \
        ${covar_info_file} \
        ${map_file} \
        ${chrom} \
        ${peak_pos} \
        ${start_pos} \
        ${end_pos} \
        ${scan1out_file}
    """

}