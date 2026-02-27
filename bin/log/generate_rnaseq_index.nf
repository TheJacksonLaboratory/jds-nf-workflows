import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){

if (params.mgi && params.annotation_source != "MGI") {
    error "The --mgi parameter is only valid when --annotation_source is set to MGI. Please set --annotation_source to MGI, or remove the --mgi parameter."
}

if (params.mgi && params.ref_table == null && params.ref_table == '') {
    error "The -mgi parameter requires that --ref_table is set. Please set --ref_table to the path of the MGI reference table, or remove the --mgi parameter."
}

def message = ""

message = """
GENERATE RNASEQ INDEX PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
-w                              ${workDir}
-c                              ${params.config}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--annotation_source             ${params.annotation_source}
--ref_fa                        ${params.ref_fa}
--ref_gtf                       ${params.ref_gtf}
--ref_gff                       ${params.ref_gff}
--ref_gff3                      ${params.ref_gff3}
--mgi                           ${params.mgi}
--ref_table                     ${params.ref_table}
--custom_gene_fasta             ${params.custom_gene_fasta}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

log.info(message)
return(message)

}
