import Logos

logo = new Logo()
println '\n'
println logo.show()


if (!(params.gen_org in ['human', 'mouse', 'other'])) {
  log.error "This workflow supports only --gen_org mouse, --gen_org human, or --gen_org other. Please retry with correct option."
  System.exit(1)
}

if (!params.fa_cds) {
  log.error "The --fa_cds parameter is required. Please provide an absolution path, or ftp download link to all transcript coding sequences."
  System.exit(1)
}

if (!(params.library_strategy in ['PE', 'SE', 'BOTH'])) {
  log.error "This workflow supports only --library_strategy PE, --library_strategy SE, or --library_strategy BOTH. Please retry with correct option."
  System.exit(1)
}

if (params.read_length != 100 && params.read_length != 75 && !params.quality_reference && params.simulate_sequencing_error) {
  log.error "If --simulate_sequencing_error is TRUE and --quality_reference is not provided, the output read length must be either 100-bp or 75-bp. Please retry with correct option."
  System.exit(1)
}

def param_log(){

def message = ""

message = """
GENERATE RNASEQ SIMREADS PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
-w                              ${workDir}
-c                              ${params.config}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--fa_cds                        ${params.fa_cds}

--library_size                  ${params.library_size}
--library_strategy              ${params.library_strategy}
--read_length                   ${params.read_length}
--fragment_length_min           ${params.fragment_length_min}
--fragment_length_max           ${params.fragment_length_max}
--fragment_length_mean          ${params.fragment_length_mean}
--fragment_length_sd            ${params.fragment_length_sd}
--simulate_sequencing_error     ${params.simulate_sequencing_error}
--quality_reference             ${params.quality_reference}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

log.info(message)
return(message)

}
