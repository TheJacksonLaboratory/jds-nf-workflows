 import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){

    if (params.gen_org != "human" && params.gen_org != "mouse") {
        error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'mouse', or 'human'" 
    }

    if (params.gen_org == 'other' && params.run_sv) {
        error "'--gen_org': When '--gen_org other' is specified, '--run_sv' must be false. SV calling is only supported for mouse and human samples."
    }

    def baseParams = """
    --workflow                      ${params.workflow}
    --gen_org                       ${params.gen_org}
    --genome_build                  ${params.genome_build}
    --gen_ver                       ${params.gen_ver}
    --csv_input                     ${params.csv_input}
    -w                              ${workDir}
    -c                              ${params.config}
    --pubdir                        ${params.pubdir}
    --ref_fa                        ${params.ref_fa}
    --ref_fa_indices                ${params.ref_fa_indices}
    --dbSNP                         ${params.dbSNP}
    --snpEff_config                 ${params.snpEff_config}
    """

    def humanParams = """
    --gold_std_indels               ${params.gold_std_indels}
    --phase1_1000G                  ${params.phase1_1000G}
    --dbNSFP                        ${params.dbNSFP}
    --cosmic                        ${params.cosmic}
    --snpEff_config                 ${params.snpEff_config}
    """

    def svParams = """
    --ref_fa_dict                   ${params.ref_fa_dict}
    --smoove_support                ${params.smoove_support}
    --exclude_regions               ${params.exclude_regions}
    --delly_exclusion               ${params.delly_exclusion}
    --delly_mappability             ${params.delly_mappability}
    --cnv_window                    ${params.cnv_window}
    --cnv_min_size                  ${params.cnv_min_size}
    --callRegions                   ${params.callRegions}
    --combined_reference_set        ${params.combined_reference_set}
    --cytoband                      ${params.cytoband}
    --sv_slop                       ${params.sv_slop}
    --sizemargin                    ${params.sizemargin}
    --min_sv_length                 ${params.min_sv_length}
    --cnv_distance_limit            ${params.cnv_distance_limit}
    --ensemblUniqueBed              ${params.ensemblUniqueBed}
    """

    def svHumanParams = """
    --dgv                           ${params.dgv}
    --thousandG                     ${params.thousandG}
    --cosmicUniqueBed               ${params.cosmicUniqueBed}
    --gap                           ${params.gap}
    --dgvBedpe                      ${params.dgvBedpe}
    --thousandGVcf                  ${params.thousandGVcf}
    --svPon                         ${params.svPon}
    --cosmicBedPe                   ${params.cosmicBedPe}
    """

    def svMouseParams = """
    --known_del                     ${params.known_del}
    --known_ins                     ${params.known_ins}
    --known_inv                     ${params.known_inv}
    """

    def logHeader = """
    WGS_SV_BAM PARAMETER LOG

    --comment: ${params.comment}

    Results Published to: ${params.pubdir}
    ________________________________________________________________________________________
    """

    def projectInfo = """
    Project Directory: ${projectDir}

    Command line call: 
    ${workflow.commandLine}
    ________________________________________________________________________________________
    """

    def msg = logHeader + baseParams + svParams
    
    if (params.gen_org == 'human') {
        msg += svHumanParams
        msg += humanParams
    } else if (params.gen_org == 'mouse') {
        msg += svMouseParams
    }

    msg += projectInfo

    log.info msg

    return(msg)
}
