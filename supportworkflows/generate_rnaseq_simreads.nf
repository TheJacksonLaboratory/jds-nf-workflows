#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "../bin/help/generate_rnaseq_simreads"
include {param_log} from "../bin/log/generate_rnaseq_simreads"
include {final_run_report} from "../bin/shared/final_run_report.nf"
include {REFORMAT_FASTA} from "../modules/bbmap/bbmap_reformat_fasta"
include {GENERATE_SIMULATED_RNA_DATA} from "../modules/r/generate_simulated_RNA_data"


workflow GENERATE_RNASEQ_SIMREADS {
    
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

    REFORMAT_FASTA(params.fa_cds)

    library_size = params.library_size instanceof String 
        ? channel.from(params.library_size.split(',').collect{it.trim().toInteger()})
        : channel.from(params.library_size.collect{it.toInteger()})

    sim_input = params.library_strategy == 'BOTH' ? REFORMAT_FASTA.out.wrap_fasta.combine(library_size).combine(channel.from(['PE', 'SE'])) : REFORMAT_FASTA.out.wrap_fasta.combine(library_size).combine(channel.from([params.library_strategy]))

    GENERATE_SIMULATED_RNA_DATA(sim_input)

}

def checkFileExists(filePath, name) {
    if (filePath && !file(filePath).exists()) {
        log.error "File not found: ${filePath} (${name})"
        exit 1
    }
}
