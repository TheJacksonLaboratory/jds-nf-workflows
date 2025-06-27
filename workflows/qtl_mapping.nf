#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/qtl_mapping.nf"
include {param_log} from "${projectDir}/bin/log/qtl_mapping.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_csv_qtl.nf"
include {DATA_QC} from "${projectDir}/modules/qtl2/data_qc.nf"
include {MAP_QTL} from "${projectDir}/modules/qtl2/map_qtl.nf"
include {RUN_PERMS} from "${projectDir}/modules/qtl2/run_perms.nf"
include {HARVEST_QTL} from "${projectDir}/modules/qtl2/harvest_qtl.nf"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

workflow QTL_MAPPING {
    
    project_ch = extract_csv(params.csv)
    project_ch.map{it -> [  it[0], 
                            it[1].covar_file, it[1].cross_file, it[1].genoprobs_file, it[1].alleleprobs_file, it[1].kinship_file, 
                            it[1].pheno_file, it[1].transformation_file] }set{data_qc_ch}

    // Data quality control
    DATA_QC(data_qc_ch)
    prob_files = DATA_QC.out.probs_files
    
    // Extract phenotype from file name 
    pheno_ch = DATA_QC.out.pheno_files.flatten().map {  file_path -> def name = file_path.name.replaceFirst(/_pheno\.csv$/, '')
                                                        [file_path, name] }

    // Combine phenotype files across prob files to do QTL mapping in parallel
    map_perm_ch = prob_files.combine(pheno_ch)

    // Map QTL
    MAP_QTL(map_perm_ch)

    // Run permutations
    RUN_PERMS(map_perm_ch)

    // Collect results
    perm_ch =   RUN_PERMS.out.perm_files.groupTuple()
    map_ch  =   MAP_QTL.out.scan1_files.groupTuple()
    harvest_ch = perm_ch.join(map_ch, by: 0)
    HARVEST_QTL(harvest_ch)

    
}