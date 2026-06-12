#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/generate_rnaseq_simreads"
include {param_log} from "${projectDir}/bin/log/generate_rnaseq_simreads"
include {final_run_report} from "${projectDir}/bin/shared/final_run_report.nf"
include {REFORMAT_FASTA} from "${projectDir}/modules/bbmap/bbmap_reformat_fasta"
include {GENERATE_SIMULATED_RNA_DATA} from "${projectDir}/modules/r/generate_simulated_RNA_data"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
message = param_log()


def checkFileExists(filePath, name) {
    if (filePath && !file(filePath).exists()) {
        log.error "File not found: ${filePath} (${name})"
        exit 1
    }
}

workflow GENERATE_RNASEQ_SIMREADS {

    REFORMAT_FASTA(params.fa_cds)

    library_size = params.library_size instanceof String 
        ? Channel.from(params.library_size.split(',').collect{it.trim().toInteger()})
        : Channel.from(params.library_size.collect{it.toInteger()})

    sim_input = params.library_strategy == 'BOTH' ? REFORMAT_FASTA.out.wrap_fasta.combine(library_size).combine(Channel.from(['PE', 'SE'])) : REFORMAT_FASTA.out.wrap_fasta.combine(library_size).combine(Channel.from([params.library_strategy]))

    GENERATE_SIMULATED_RNA_DATA(sim_input)

}
