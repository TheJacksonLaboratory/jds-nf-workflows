import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){
log.info """
GENERATE RNASEQ INDEX PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
-w                              ${workDir}
-c                              ${params.config}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--ref_fa                        ${params.ref_fa}
--ref_gtf                       ${params.ref_gtf}
--ref_gff                       ${params.ref_gff}
--ref_gff3                      ${params.ref_gff3}
--custom_gene_fasta             ${params.custom_gene_fasta}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
}
