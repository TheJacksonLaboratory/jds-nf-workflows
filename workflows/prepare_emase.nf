#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "../bin/help/prepare_emase.nf"
include {param_log} from "../bin/log/prepare_emase.nf"
include {final_run_report} from "../bin/shared/final_run_report.nf"
include {EMASE_PREPARE_EMASE} from "../modules/emase/emase_prepare_emase"
include {BOWTIE_BUILD} from "../modules/bowtie/bowtie_build"
include {CLEAN_TRANSCRIPT_LISTS} from "../modules/python/clean_prepEmase_transcriptList"

// main workflow
workflow PREPARE_EMASE {

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

    // Prepare emase reference, given list of genomes and gtf files. 
    EMASE_PREPARE_EMASE()
    BOWTIE_BUILD(EMASE_PREPARE_EMASE.out.pooled_transcript_fasta, 'bowtie.transcripts')
    // clean transcript lists to add transcripts absent from certain haplotypes.
    CLEAN_TRANSCRIPT_LISTS(EMASE_PREPARE_EMASE.out.pooled_transcript_info)
}
