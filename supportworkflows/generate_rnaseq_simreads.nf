#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/generate_rnaseq_simreads"
include {param_log} from "${projectDir}/bin/log/generate_rnaseq_simreads"
include {REFORMAT_FASTA} from "${projectDir}/modules/bbmap/bbmap_reformat_fasta"
include {GENERATE_SIMULATED_RNA_DATA} from "${projectDir}/modules/r/generate_simulated_RNA_data"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
message = param_log()

// Save params to a file for record-keeping
workflow.onComplete {
    final_run_report(message)
}

def checkFileExists(filePath, name) {
    if (filePath && !file(filePath).exists()) {
        log.error "File not found: ${filePath} (${name})"
        exit 1
    }
}


workflow GENERATE_RNASEQ_SIMREADS {

    REFORMAT_FASTA(params.fa_cds)
    GENERATE_SIMULATED_RNA_DATA(REFORMAT_FASTA.out.wrap_fasta)

}

