import Logos

logo = new Logo()
println '\n'
println logo.show()

if (!params.fasta) {
  log.error "The --fasta parameter is required. Please provide the path to reference fasta."
  System.exit(1)
}

if (!params.read_type || !(params.read_type in ['PE', 'SE'])) {
  log.error "The --read_type parameter is required. Please specify 'PE' for paired-end or 'SE' for single-end reads."
  System.exit(1)
}

def param_log(){

def message = ""

message = """
GENERATE WGS SIMREADS PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
-w                              ${workDir}
-c                              ${params.config}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--sampleID                      ${params.sampleID}
--fasta                         ${params.fasta}
--chrom_list                    ${params.chrom_list}
--coverage                      ${params.coverage}
--read_length                   ${params.read_length}
--error_rate                    ${params.error_rate}
--read_type                     ${params.read_type}
${params.read_type == 'PE' ? "--insert_size                   ${params.insert_size}" : ""}
${params.read_type == 'PE' ? "--insert_size_sd                ${params.insert_size_sd}" : ""}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

log.info(message)
return(message)

}
