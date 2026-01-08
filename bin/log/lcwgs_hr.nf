import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){

if (params.gen_org == "human" || params.gen_org == "other") {
        error "'--gen_org': \"${params.gen_org}\" is not valid, only supported option is 'mouse' for lcwgs_hr workflow." 
}

if (params.covar_file == null) {
        error "'--covar_file': \"${params.covar_file}\" is not valid for lcwgs_hr workflow." 
}

    def baseParams = """
    --workflow                      ${params.workflow}
    --gen_org                       ${params.gen_org}
    --genome_build                  ${params.genome_build}
    --read_type                     ${params.read_type}
    --sample_folder                 ${params.sample_folder}
    --pattern                       ${params.pattern}
    --extension                     ${params.extension}
    --concat_lanes                  ${params.concat_lanes}
    --csv_input                     ${params.csv_input}
    --download_data                 ${params.download_data}
    --library_type                  ${params.library_type}
    --gridfile                      ${params.gridfile}
    --covar_file                    ${params.covar_file}
    --cross_type                    ${params.cross_type}
    --smooth_window                 ${params.smooth_window}
    -w                              ${workDir}
    -c                              ${params.config}
    --pubdir                        ${params.pubdir}
    --ref_fa                        ${params.ref_fa}
    --ref_fa_indices                ${params.ref_fa_indices}
    --ref_haps_dir                  ${params.ref_haps_dir}
    --unqualified_perc              ${params.unqualified_perc}
    --detect_adapter_for_pe         ${params.detect_adapter_for_pe}
    --mismatch_penalty              ${params.mismatch_penalty}
    """


    def logHeader = """
    LCWGS_HR PARAMETER LOG

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

    def msg = logHeader + baseParams
    
    msg += projectInfo

    log.info msg
}
