import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){

def message = ""

message = """
HAPLOTYPE RECONSTRUCTION PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--csv_input                     ${params.csv_input}
--n_perms                       ${params.n_perms}
--rerun                         ${params.rerun}
--correct_ids                   ${params.correct_ids}
--remove_markers                ${params.remove_markers}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--gigamuga_reference            ${params.gigamuga_reference}
--gm_gmaps                      ${params.gm_gmaps}
--gm_pmaps                      ${params.gm_pmaps}
--gm_cc_do_founder_genotypes    ${params.gm_cc_do_founder_genotypes}
--gm_cc_do_allele_codes         ${params.gm_cc_do_allele_codes}
--gm_ann = ${params.gm_ann}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

log.info(message)
return(message)

}
