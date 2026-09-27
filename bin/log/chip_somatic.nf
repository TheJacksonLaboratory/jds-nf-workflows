def param_log() {

if (params.gen_org != "mouse" && params.gen_org != "human") {
    error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'human' or 'mouse'"
}

def msg = """
CHIP_SOMATIC PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________

--workflow                      ${params.workflow}
--csv_input                     ${params.csv_input}
--pubdir                        ${params.pubdir}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--read_type                     ${params.read_type}

--quality_phred                 ${params.quality_phred}
--unqualified_perc              ${params.unqualified_perc}
--detect_adapter_for_pe         ${params.detect_adapter_for_pe}
--trim_poly_g                   ${params.trim_poly_g}
--trim_poly_x                   ${params.trim_poly_x}
${params.trim_poly_x ? "--poly_x_min_len                ${params.poly_x_min_len}" : ""}

--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}

--target_gatk                   ${params.target_gatk}
--target_picard                 ${params.target_picard}
--bait_picard                   ${params.bait_picard}

--germline_resource_vcf         ${params.germline_resource_vcf}
--pon_ref                       ${params.pon_ref}
--genotype_pon                  ${params.genotype_pon}
--genotype_germline             ${params.genotype_germline}
--contam_ref                    ${params.contam_ref}
______________________________________________________
"""

def projectInfo = """
Project Directory: ${projectDir}

Command line call:
${workflow.commandLine}
______________________________________________________
"""

msg += projectInfo

log.info msg

return(msg)

}
