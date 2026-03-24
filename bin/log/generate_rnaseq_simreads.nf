import Logos

logo = new Logo()
println '\n'
println logo.show()


if (!(params.gen_org in ['human', 'mouse'])) {
  log.error "This workflow supports only --gen_org mouse or --gen_org human. Please retry with correct option."
  System.exit(1)
}

if (!params.fa_cds) {
  log.error "The --fa_cds parameter is required. Please provide ftp download link to all transcript coding sequences resulting from Ensembl genes."
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
--annotation_source             ${params.annotation_source}
--fa_cds                        ${params.fa_cds}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

log.info(message)
return(message)

}
