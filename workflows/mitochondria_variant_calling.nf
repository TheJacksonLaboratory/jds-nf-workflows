#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "../bin/help/mitochondria_variant_calling.nf"
include {param_log} from "../bin/log/mitochondria_variant_calling.nf"
include {final_run_report} from "../bin/shared/final_run_report.nf"
include {extract_csv_bam} from "../bin/shared/extract_csv_bam.nf"
include {MT_VARIANT_CALLING} from "../subworkflows/mt_variant_calling"

// main workflow
workflow MITOCHONDRIA_VARIANT_CALLING {
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

    bam_ch = extract_csv_bam(file(params.csv_input, checkIfExists: true))
    input_ch = bam_ch.map{it -> [it[0], file(it[2]), file(it[3])]}

    MT_VARIANT_CALLING(input_ch)
    // workflow found in: subworkflows/mt_variant_calling.nf
    // workflow run as subworkflow due to re-use in WGS and PTA (possible) workflows. 
}
