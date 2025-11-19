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
include {QTL_EFFECTS} from "${projectDir}/modules/qtl2/qtl_effects.nf"
include {SUMMARIZE_QTL_EFFECTS} from "${projectDir}/modules/qtl2/summarize_qtl_effects.nf"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

workflow QTL_MAPPING {
    
    project_ch = extract_csv(params.csv_input)
    project_ch.map{it -> [  it[0],
                            it[1].covar_file, it[1].map_file, it[1].genoprobs_file, it[1].alleleprobs_file, it[1].kinship_file, 
                            it[1].pheno_file, it[1].covar_info_file] }set{data_qc_ch}
    data_qc_ch.map{it -> [  it[0], it[2] ] }.set{map_file_ch}
    
    // Data quality control
    DATA_QC(data_qc_ch.map{it -> [  it[0], it[1], it[3], it[4], it[5], it[6], it[7] ] })
    prob_files = DATA_QC.out.probs_files
    
    // Extract phenotype from file name
    pheno_ch        = DATA_QC.out.pheno_files.flatten().map {  file_path -> def name = file_path.name.replaceFirst(/_pheno\.csv$/, '')
                                                        [file_path, name] }
    covar_info_ch   = DATA_QC.out.covar_info_files.flatten().map {  file_path -> def name = file_path.name.replaceFirst(/_covar_info\.csv$/, '')
                                                        [file_path, name] }
    
    pheno_covar_info_ch = pheno_ch.join(covar_info_ch, by: 1)
    
    // Combine phenotype files and covar info files across prob files to do QTL mapping in parallel
    map_perm_ch = prob_files.combine(pheno_covar_info_ch)
    map_perm_ch = map_perm_ch.combine(map_file_ch, by: 0)
    
    // Map QTL
    MAP_QTL(map_perm_ch)

    // Run permutations
    RUN_PERMS(map_perm_ch)

    // // Collect results
    perm_ch =   RUN_PERMS.out.perm_files.groupTuple()
    map_ch  =   MAP_QTL.out.scan1_files.groupTuple()
    harvest_ch = perm_ch.join(map_ch, by: 0)
    HARVEST_QTL(harvest_ch)

    // Make channels for each QTL
    HARVEST_QTL.out.qtl_table.splitCsv(header: true).map{
        it -> [ [it[0], it[1].lodcolumn], [it[1].chr, it[1].pos, it[1].ci_lo, it[1].ci_hi] ]
    }.set{peaks_ch}

    // Collect scan1 files for effect estimation
    MAP_QTL.out.scan1_files
               .map{
                    it -> [ [it[0], it[1]], it[2] ]
               }.set{scan1_ch}

    // Join peaks with probs files for effect estimation
    MAP_QTL.out.probs_files
               .map{
                    it -> [ [it[0], it[5]], [it[1], it[2], it[3], it[4], it[6], it[7], it[8]] ] // sets the project id and phenotype as combine key
               }.combine(peaks_ch, by: 0)
               .combine(scan1_ch, by: 0)
               .map{
                    it -> [it[0][0], it[0][1], // project id and phenotype
                           it[1][0], it[1][1], it[1][2], it[1][3], it[1][4], it[1][5], it[1][6], // prob files
                           it[2][0], it[2][1], it[2][2], it[2][3], // peak info
                           it[3]] // scan1 file
               }.set{probs_peaks_ch}

    // Estimate QTL effects
    QTL_EFFECTS(probs_peaks_ch)
    
    // Summarize QTL effects
    SUMMARIZE_QTL_EFFECTS(QTL_EFFECTS.out.qtl_peaks_files.groupTuple())
}